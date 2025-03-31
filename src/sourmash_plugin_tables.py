"""\
The tables plugin for sourmash output files

    Reads and processes any number of sourmash output files to create a table of values for `match_name` per `query_name`.

    Parameters:
        filename (str): Path to the sourmash file (gather or prefetch CSV file).
        taxononmy_file (str): Path to the taxonomy file containing the lineage ranks for the sourmash files.
        lineage_rank (str): The column header that corresponds to the lineage rank for compressing the sourmash files.
        column (str): The column header of the sourmash files to populate the values of the table.
        output (str): The output file name.
        gzip: Use for compressed file type.
        output_format (str): Either "dense" or "sparse" to specify the output format.

    Returns:
        pl.DataFrame converted to CSV file:
            - output_format:
                - Dense: A DataFrame with `query_name` as rows and `match_name` as columns.
                - Sparse: A DataFrame with rows indicating `query_name`, `match_name`, and their value.
            - taxonomy_file and lineage_rank:
                - A summation of the values for all lower order lineages that populate that row.
"""

usage="""
  sourmash scripts gather_tables gather-dir/*.gather.csv --output gather.csv
  sourmash scripts prefetch_tables  prefetch-dir/*.prefetch.csv --output prefetch.csv
"""

epilog="""
See https://github.com/sourmash-bio/sourmash_plugin_tables for more examples.

Need help? Have questions? Ask at http://github.com/sourmash-bio/sourmash/issues!
"""

import argparse
import sourmash

from sourmash.index import LinearIndex
from sourmash.logging import debug_literal
from sourmash.plugins import CommandLinePlugin

from sourmash.save_load import (Base_SaveSignaturesToLocation,
                                _get_signatures_from_rust)

import sys
import polars as pl
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

class Command_Prefetch_Tables(CommandLinePlugin):
    command = 'prefetch_tables'             # 'scripts <command>'
    description = __doc__       # output with -h
    usage = usage               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, parser_prefetch):
        super().__init__(parser_prefetch)

        # Subparser for 'prefetch'
        parser_prefetch.add_argument('filenames', nargs='+', help="List of input sourmash prefetch files to combine.")
        parser_prefetch.add_argument('-t', '--taxonomy-file', '--taxonomy', nargs='?', metavar='FILE', default=None, help="The Sourmash taxonomy file that corresponds with the database that generated the sourmash prefetch files.")
        parser_prefetch.add_argument('-l', '--lineage-rank', '--lineage', default='s__species', help="Compress the tables by user-defined taxonomic lineage rank associated with sourmash taxonomic file")
        parser_prefetch.add_argument('-c', '--column', type=str, default='jaccard', help="The numerical column from 'prefetch' to populate the table values (default: 'jaccard').\nVisit https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#id23 for more information.")
        parser_prefetch.add_argument('-o', '--output', required=True, help="Path to save the combined output CSV file.")
        parser_prefetch.add_argument('-f', '--format', choices=['dense', 'sparse'], default='dense', help="Output file structure: dense or sparse OTU.")
        parser_prefetch.add_argument('-z', '--gzip', action='store_true', help="Compress the output file into a .gz file type.")

        debug_literal('RUNNING cmd_prefetch_tables.__init__')

    def main(self, args):
        super().main(args)

        tables_main(args)

        print('RUNNING cmd', self, args)

class Command_Gather_Tables(CommandLinePlugin):
    command = 'gather_tables'             # 'scripts <command>'
    description = __doc__       # output with -h
    usage = usage               # output with no args/bad args as well as -h
    epilog = epilog             # output with -h
    formatter_class = argparse.RawTextHelpFormatter # do not reformat multiline

    def __init__(self, parser_gather):
        super().__init__(parser_gather)

        # Subparser for 'gather'
        parser_gather.add_argument('filenames', nargs='+', help="List of input sourmash gather files to combine.")
        parser_gather.add_argument('-t', '--taxonomy-file', '--taxonomy', nargs='?', metavar='FILE', default=None, help="The Sourmash taxonomy file that corresponds with the database that generated the sourmash gather files.")
        parser_gather.add_argument('-l', '--lineage-rank', '--lineage', default='s__species', help="Compress the tables by user-defined taxonomic lineage rank associated with sourmash taxonomic file")
        parser_gather.add_argument('-c', '--column', type=str, default='f_unique_weighted', help="The numerical column from 'gather' to poplate the table values (default: 'f_unique_weighted').\nVisit https://sourmash.readthedocs.io/en/latest/classifying-signatures.html#id22 for more information.")
        parser_gather.add_argument('-o', '--output', required=True, help="Path to save the combined output CSV file.")
        parser_gather.add_argument('-f', '--format', choices=['dense', 'sparse'], default='dense', help="Output file structure: dense or sparse OTU.")
        parser_gather.add_argument('-z', '--gzip', action='store_true', help="Compress the output file into a .gz file type.")

        debug_literal('RUNNING cmd_gather_tables.__init__')

    def main(self, args):
        super().main(args)

        tables_main(args)

        print('RUNNING cmd', self, args)


def process_file(filename, taxdb, output_format="dense", column_selection='intersect_bp', lineage_rank='species'):
    """
    Reads and processes any number of sourmash files to create a matrix of `match_name` values per `query_name`.

    Parameters:
        filename (str): Path to the sourmash file.
        output_format (str): Either "dense" or "sparse" to specify the output format.

    Returns:
        pl.DataFrame: 
            - Dense: A DataFrame with `query_name` as rows and `match_name` as columns.
            - Sparse: A DataFrame with rows indicating `query_name`, `match_name`, and their value.

            Note: column names are the basename of the file prior to a '.' I.e. ERR1234567.sig will populate the ERR1234567 column
    """

    try:
        df = pl.read_csv(
            filename,
            separator=',',
            has_header=True
        )

        if 'name' in df.columns:
            df = df.rename({'name': 'match_name'})

        # check for and select the columns of interest for the table
        required_columns = [column_selection, 'query_name', 'match_name']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"Missing required columns in file: {filename}")

        column_names = ['query_name', 'match_name', column_selection]
        df = df.select(column_names)

        if taxdb is not None:
            # find match_name that is in tax 
            df = df.with_columns(
                    pl.col('match_name').str.split(' ').list.get(0).alias('match_ident')
                    )
            merged_df = df.join(taxdb, left_on='match_ident', right_on='ident', how='left')
            missing_df = merged_df.filter(pl.col(lineage_rank).is_null()).drop(['match_ident', lineage_rank])

            column_names.append(lineage_rank)

            tax_df = (merged_df
                    .select(column_names)
                    .rename({lineage_rank: f'match_name_{lineage_rank}'})
                    .drop_nulls()
                    .group_by([f'match_name_{lineage_rank}', 'query_name'])
                    .agg(pl.sum(column_selection))
                    )
            df = tax_df


        #print(df, merged_df, missing_df, tax_df)

        # create the table for a dense or sparse format
        if output_format == "dense" :
            dense_matrix = df.pivot(
                values = column_selection,
                index = f"match_name_{lineage_rank}" if taxdb is not None and f"match_name_{lineage_rank}" in df.columns else "match_name",
                columns = "query_name",
            ).fill_null(0)  # replace null values with 0
            return dense_matrix

        elif output_format == "sparse":
            return df

        else:
            raise ValueError(f"Invalid output_format: {output_format}. Use 'dense' or 'sparse'.")

    # What specific exceptions should I expect?
    except Exception as e:
        raise RuntimeError(f"Error processing file {filename}: {e}")

def process_file_with_format(args):
    """
    A wrapper for process_file function to allow multiprocessing with additional arguments.
    """
    filename, taxdb, output_format, column_selection, rank = args
    return process_file(filename, taxdb, output_format=output_format, column_selection=column_selection, lineage_rank=rank)

def tables_main(args):

    args

    if args.taxonomy_file:
        print(f"loading taxonomies from {args.taxonomy_file}")
        #taxdb = sourmash.tax.tax_utils.MultiLineageDB.load([args.taxonomy_file])
        taxdb = pl.read_csv(args.taxonomy_file, separator=',', has_header=True)
        print(f"found {len(taxdb)} identifiers in taxdb.")
        print(taxdb)
        #ident_set = set(taxdb["ident"].to_list())
        #print(len(ident_set))

        # Parallel processing all files (list a tuple of filename with output structure, send that list to process_file_with_format which runs the process_file function)
        file_format_args = [(filename, taxdb, args.format, args.column, args.lineage_rank) for filename in args.filenames]

        with ProcessPoolExecutor() as executor:
            dfs = list(executor.map(process_file_with_format, file_format_args))

    else:
        # Parallel processing all files (list a tuple of filename with output structure, send that list to process_file_with_format which runs the process_file function)
        file_format_args = [(filename, args.taxonomy_file, args.format, args.column, args.lineage_rank) for filename in args.filenames]

        with ProcessPoolExecutor() as executor:
            dfs = list(executor.map(process_file_with_format, file_format_args))
    print(dfs)

    # Combine all DataFrames
    if args.format == "dense":
        # Outer join to create a dense matrix with all 'name' columns
        # https://docs.pola.rs/api/python/dev/reference/api/polars.concat.html#polars.concat
        combined_df = pl.concat(dfs, how='align')

        combined_df = combined_df.fill_null(0)
        print(combined_df)
        # Create a list of data and numeric columns
        #data_cols = [col for col in combined_df.columns if col not in {'match_name', f'match_name_{args.lineage_rank}'}]
        numeric_cols = combined_df.select(pl.col(pl.Int64) | pl.col(pl.Float64) | pl.col(pl.UInt64)).columns

#        # Compute the maximum percentage across all input files for each taxonomic name
#        combined_df = combined_df.with_columns(
#            pl.concat_list(pl.col(numeric_cols)).arr.max().alias('max')
#        )
#
#        # Count values > 0 in numeric columns i.e what is the total count of genome across samples
#        combined_df = combined_df.with_columns(
#            pl.sum_horizontal((pl.col(numeric_cols) > 0).cast(pl.Int64)).alias("count")
#        )
#
#        # Add 'total' column: The sum of the selected columns
#        combined_df = combined_df.with_columns(
#            pl.concat_list(pl.col(numeric_cols)).arr.sum().alias('total')
#        )

    else:  # Sparse format
        # Concatenate DataFrames without merging on 'name' to maintain sparse format
        combined_df = pl.concat(dfs)
    
        # Optionally sort by the percentage column ('percent' or similar)
        if 'f_unique_weighted' in combined_df.columns:
            combined_df = combined_df.sort('f_unique_weighted', descending=True)

    # What should be the final output upon completion?
    print(combined_df)

    print(combined_df.shape)
    print(combined_df.schema)
    #print(combined_df.describe())
    print(combined_df.sample(fraction=0.01))

    if args.gzip:
        output = args.output + '.gz'

        with gzip.open(output, 'wt', encoding='UTF-8') as fp:
            combined_df.write_csv(fp)
       # # Test datasets
       # combined_df = pl.DataFrame({
       #     "match_name": ["apple", "banana", "orange", "grape", "melon"]
       # })
       # combined_df = pl.DataFrame({
       #     "match_name": ["apple,banana", "orange", "grape,kiwi", "melon"]
       # })

        # Print something to help look at the output in the terminal
        if combined_df["match_name"].str.contains(",").any() or combined_df[f"match_name_{args.lineage_rank}"].str.contains(",").any():
            print(f"""Consider running `gzip -cd {output} | sed -E 's/"([^"]*),([^"]*)"/"\\1|\\2"/g' | column -s, -t | less -S` to see the full table.""")
        else:
            print(f"Consider running `gzip -cd {output} | column -s, -t | less -S` to see full table.")

    else:
        combined_df.write_csv(args.output)

        if "match_name" in combined_df.columns:
            if combined_df["match_name"].str.contains(",").any():
                print(f"""Consider running `sed -E 's/"([^"]*),([^"]*)"/"\\1|\\2"/g' {args.output} | column -s, -t | less -S` to see the full table.""")
            else:
                print(f"Consider running `cat {args.output} | column -s, -t | less -S` to see full table.")

        elif f"match_name_{args.lineage_rank}" in combined_df.columns:
            if combined_df[f"match_name_{args.lineage_rank}"].str.contains(",").any():
                print(f"""Consider running `sed -E 's/"([^"]*),([^"]*)"/"\\1|\\2"/g' {args.output} | column -s, -t | less -S` to see the full table.""")
            else:
                print(f"Consider running `cat {args.output} | column -s, -t | less -S` to see full table.")

