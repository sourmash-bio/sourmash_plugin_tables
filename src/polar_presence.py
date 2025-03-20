#! /usr/bin/env python

import sys
from itertools import islice
import argparse
import sourmash
import polars as pl

import sourmash_utils

def read_file_and_separate(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    header = lines[0].strip()
    idx = [line.strip() for line in lines[1:]]

    return header, idx


def main():
    p = argparse.ArgumentParser()

    p.add_argument('ranktable', help="Input csv containing classified hashes for specified organism")

    # Required input: Sketches file (containing MinHash sketches for metagenomes)
    p.add_argument('sketches', help="Input file with sketches to process")

    # Required output: Output file to store hash presence/absence info
    p.add_argument('-o', '--output', required=True, help="Output CSV/Parquet file")

    # Optional: File to filter which samples to include based on names
    p.add_argument('--filter-samples', default=None, help="Optional file with sample names to include")

    # Choose output format: CSV or Parquet
    p.add_argument('--format', choices=['csv', 'parquet'], default='csv',
                   help="Output file format: 'csv' or 'parquet' (default: csv)")

    p.add_argument('-c','--collapse-columns', nargs="*", help='Collapse the polars dataframe by the header of each text file')

    sourmash_utils.add_standard_minhash_args(p)

    args = p.parse_args()

    print(f"Loading ranktable CSV '{args.ranktable}'...")
    ranktable_df = pl.read_csv(args.ranktable)
    print(ranktable_df)
    print(sys.getsizeof(ranktable_df), 'bytes')

#    classify_d = read_ranktable_csv(args.ranktable)
#    print(f"Loaded {len(classify_d)} hashvals... downsampling soon.")
#    print(sys.getsizeof(classify_d), 'bytes')
#    print(dict(islice(classify_d.items(), 0, 5)))
#    print('     ...    ')
#    last_five = dict(islice(reversed(classify_d.items()), 5))
#    last_five = dict(reversed(list(last_five.items())))
#    print(last_five)

    hashvals_l = ranktable_df['hashval'] # extract as a series or polars list
    print(f"Loaded {len(hashvals_l)} hashvals... downsampling soon.")
    print(hashvals_l)
    print(sys.getsizeof(hashvals_l), 'bytes')

    select_mh = sourmash_utils.create_minhash_from_args(args)
    print(f"selecting sketches: {select_mh}")

    # Load the samples
    print(f"loading sketches from file '{args.sketches}'")
    idx = sourmash_utils.load_index_and_select(args.sketches, select_mh)

    print(f"found {len(idx)} metagenomes")

    print(idx)

    # Create an empty MinHash to hold hashes for query
    query_minhash = next(iter(idx.signatures())).minhash.copy_and_clear()

    print(query_minhash)
    # Add all hashes from ranktable_df to the query MinHash
    query_minhash.add_many(hashvals_l.to_list())
    print(query_minhash)
    # Downsample query MinHash if `scaled` is provided
    if args.scaled:
        query_minhash = query_minhash.downsample(scaled=args.scaled)
    else:
        args.scaled = query_minhash.scaled


    # Optional: Read list of sample names to filter
    filter_by_name = None
    if args.filter_samples:
        filter_by_name = set([x.strip() for x in open(args.filter_samples)])

    # Prepare to store presence/absence data
    #presence_data = []
    presence_df = pl.DataFrame({"hashval": hashvals_l})


    # Track number of skipped sketches
    n_skipped = 0

    # Iterate through all sketches to check hash presence
    for n, metag_ss in enumerate(idx.signatures()):
        metag_name = metag_ss.name

        # Print progress every 10 sketches
        if n and n % 10 == 0:
            print('...', n, 'and', sys.getsizeof(presence_df), 'bytes')

        # Skip sketches if not in the filter list
        if filter_by_name and metag_name not in filter_by_name:
            n_skipped += 1
            continue

        # Downsample the current sketch to the specified scale
        metag_mh = metag_ss.minhash.downsample(scaled=args.scaled)

        # Get present hashes in the current sketch
        #metag_hashes = set(metag_mh.hashes)
        metag_hashes = pl.Series(list(metag_mh.hashes))


        #print(f"Checking hash presence for {metag_name}")
        #print(f"Number of hashes in metag_hashes: {len(metag_hashes)}")
        #print(f"First 10 hashvals_l: {hashvals_l.head(10).to_list()}")
        #print(f"First 10 metag_hashes: {metag_hashes.head(10).to_list()}")

        # Create a presence/absence list for the current sample
#        presence_data.append(
#            pl.DataFrame({
#                "hashval": hashvals_l,
#                metag_name: hashvals_l.is_in(metag_hashes).cast(pl.Int32)
#            })
#        )
        # Create presence/absence column for the current sample
        presence_column = hashvals_l.is_in(metag_hashes).cast(pl.Int32).alias(metag_name)

        # Append the new column to the existing presence_df
        presence_df = presence_df.with_columns(presence_column)
        #print(presence_df)
    print(presence_df.shape)
    #print(presence_data)
    # Merge all presence/absence dataframes
    #presence_df = pl.concat(presence_data, how="align")

    if args.collapse_columns:
        collapse_df = pl.DataFrame({"hashval": hashvals_l})
        print(f"Processing the following files: {args.collapse_columns}")
        for file in args.collapse_columns:
            print(f"Processing: {file}")
            h, i = read_file_and_separate(file)
            print(f"First Line: {h}")
            print(f"Remaining Lines: {i}")

            # Check that all `i` columns exist in presence_df
            existing_columns = [col for col in i if col in presence_df.columns]
            print(existing_columns) 
#            # Sum the values across the selected columns (ignoring `hashval`)
#            new_column = pl.col(existing_columns).sum().alias(f"{h}")
#            
#            # Add the new summed column to the dataframe
#            collapse_df = presence_df.with_columns(new_column)
#    
#            print(f"Updated DataFrame with {h}:")
#            print(collapse_df)
 
            # Only proceed if matching columns exist
            if existing_columns:
                # Sum values of matching columns
                new_column = (
                    presence_df.select(
                        sum(pl.col(col) for col in existing_columns).alias(f"{h}")
                    ).to_series()
                )

                print(new_column)   
                collapse_df = collapse_df.with_columns(new_column)

                print(f"Updated DataFrame with {h}:")
                print(collapse_df)
            else:
                print(f"No matching columns found for {file}, skipping...")
   
        # Print final dataframe after processing all files
        print("Final collapse_df:")
        print(collapse_df)

        sum_df = collapse_df.select([
            pl.col(col).sum().alias(col) for col in collapse_df.columns if col != "hashval"
        ])
        print(sum_df)
    
        print(pl.concat([pl.DataFrame({'hashval': ['count']}), sum_df], how='horizontal'))
    
        # Count values > 0 in numeric columns i.e what is the total count of genome across samples
        hor_sum_df = collapse_df.with_columns(
            pl.sum_horizontal((pl.col(pl.Int32)).cast(pl.Int64)).alias("count")
        )
    
        print(hor_sum_df)
    
        df = pl.concat([pl.DataFrame({'hashval': ['count']}), sum_df], how='horizontal')
    
        collapse_df = pl.concat([collapse_df.with_columns(pl.col("hashval").cast(pl.String)), df])
    
        hor_sum_df = collapse_df.with_columns(
            pl.sum_horizontal((pl.col(pl.Int32)).cast(pl.Int64)).alias("count")
        )
    
        print(hor_sum_df)


    else:
        print("No files provided. Nothing to process.")


    # Remove duplicate `hashval` columns after horizontal concatenation
#    presence_df = presence_df.unique(subset=["hashval"], keep="first")
    print(presence_df)

    sum_df = presence_df.select([
        pl.col(col).sum().alias(col) for col in presence_df.columns if col != "hashval"
    ])
    print(sum_df)

    print(pl.concat([pl.DataFrame({'hashval': ['count']}), sum_df], how='horizontal'))

    # Count values > 0 in numeric columns i.e what is the total count of genome across samples
    hor_sum_df = presence_df.with_columns(
        pl.sum_horizontal((pl.col(pl.Int32)).cast(pl.Int64)).alias("count")
    )

    print(hor_sum_df)

    df = pl.concat([pl.DataFrame({'hashval': ['count']}), sum_df], how='horizontal')

    presence_df = pl.concat([presence_df.with_columns(pl.col("hashval").cast(pl.String)), df])

    hor_sum_df = presence_df.with_columns(
        pl.sum_horizontal((pl.col(pl.Int32)).cast(pl.Int64)).alias("count")
    )

    print(hor_sum_df)


    # Save the presence/absence matrix to the specified format
    if args.format == 'csv':
        presence_df.write_csv(args.output)
    elif args.format == 'parquet':
        presence_df.write_parquet(args.output)

    print(f"Results written to {args.output}")
    print(f"Skipped {n_skipped} samples.")

if __name__ == '__main__':
    sys.exit(main())
