#! /usr/bin/env python

import argparse
import polars as pl
import os
import sys

def read_input_file(input_file: str) -> pl.DataFrame:
    """
    Read the input CSV or TSV file (compressed or uncompressed) into a Polars DataFrame.
    """

    if input_file.endswith('.gz'):
        # remove last 3 characters '.gz' and get file extension
        file_extension = os.path.splitext(input_file[:-3])[1].lower()
    else:
        file_extension = os.path.splitext(input_file)[1].lower()

    print("File extension: ", file_extension)

    # Determine separator based on file type (CSV or TSV)
    if file_extension == '.csv':
        separator = ',' 
    elif file_extension == '.tsv':
        separator = '\t'
    else:
        raise ValueError("Unsupported file format. Please provide a .csv or .tsv file.")

    return pl.read_csv(input_file,
                       separator=separator,
                       null_values=["NULL", "NA"],  # Add any other strings representing NULL
                       infer_schema_length=10000)  # You can increase this if necessary)

def main():
    # Set up the argument p
    p = argparse.ArgumentParser(description="Convert CSV or TSV file to Parquet format.")
    p.add_argument("input_file", help="Input CSV or TSV file (can be compressed).")
    p.add_argument("output_file", help="Output Parquet file.")
    p.add_argument("-v", "--verbose", action="store_true", help="Display max column for dataframe")

    args = p.parse_args()

    # Ensure the output file has a .parquet extension
    if not args.output_file.lower().endswith('.parquet'):
        print("Error: Output file must have a .parquet extension.")
        exit(1)

    # Read the input file into a Polars DataFrame
    print(f"Reading input file: {args.input_file}")
    df = read_input_file(args.input_file)

    if args.verbose:
        with pl.Config(tbl_cols=-1):
            print(df)
            print(df.describe())
            print(df.dtypes)
    else:
        print(df)
        print(df.describe())
        print(df.dtypes)

    # Write the DataFrame to a Parquet file
    print(f"Saving to Parquet file: {args.output_file}")
    df.write_parquet(args.output_file)
    
    print(f"File saved successfully: {args.output_file}")

if __name__ == "__main__":
    sys.exit(main())
