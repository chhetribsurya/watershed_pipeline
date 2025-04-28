#!/usr/bin/env python

import argparse
import os
import pandas as pd

def main(args):
    # Read the input file into a DataFrame
    df = pd.read_csv(args.input, sep='\t', header=None, names=['subjid', 'population'])

    # Ensure the output directory exists
    out_dir = os.path.join(args.output_dir, "outputs", "pop_ids")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Output directory: {out_dir}")

    # Iterate through each unique population and write to a separate file
    for pop in df['population'].unique():
        pop_df = df[df['population'] == pop]
        pop_file = os.path.join(out_dir, f"{pop}_ids.txt")
        print(f"Writing file: {pop_file}")
        #pop_df.to_csv(pop_file, sep='\t', header=False, index=False)
        pop_df.to_csv(f"{pop}_ids.txt", sep='\t', header=False, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split population file into separate files based on population labels.")
    parser.add_argument("--input", required=True, help="Input file with subjid and population labels.")
    parser.add_argument("--output_dir", required=True, help="Output directory for the population files.")
    args = parser.parse_args()
    main(args)
