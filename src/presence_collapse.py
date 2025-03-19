#! /user/bin/env python

import sys
import argparse
import polars as pl

def main():
    p = argparse.ArgumentParser()

    p.add_argument('presence_table', help='Input CSV or Parquet file from polars_presence.py output')
    p.add_argument('-', help='Input CSV which defines sample name by category')

    args = p.parse_args()

    
