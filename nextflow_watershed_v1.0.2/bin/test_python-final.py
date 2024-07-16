import pandas as pd
import argparse
import os

def process_file(input_file, output_file, output_dir, population):
    """
    Process the input TSV file to count non-zero rows, sum columns, drop zero-sum columns,
    and save results to specified output files.

    Args:
        input_file (str): Path to the input TSV file.
        output_file (str): Path where the processed TSV file will be saved.
        output_dir (str): Directory where the feature counts and dropped columns files will be saved.
        population (str): Population name to be used in output file names.
    """
    # Read the input file with correct handling of "NA"
    df = pd.read_csv(input_file, sep='\t', dtype={"N2pair": pd.Int64Dtype()})
    df['N2pair'] = df['N2pair'].astype(str).replace({'nan': 'NA', '<NA>': 'NA'})

    # Initialize lists to store feature counts and dropped columns
    feature_counts = []
    dropped_columns = []

    # Iterate over each column except the first two and the last column
    for column in df.columns[2:-1]:
        total_sum = df[column].sum()
        non_zero_count = (df[column] != 0).sum()
        
        # Print the total sum for the column
        print(f"Column: {column}, Total Sum: {total_sum}, Non-zero Count: {non_zero_count}")
        
        # Save the feature counts
        feature_counts.append([column, total_sum, non_zero_count])
        
        # Check if the total sum is zero and drop the column if true
        if total_sum == 0:
            dropped_columns.append(column)
            df = df.drop(columns=[column])

    # Save the feature counts to a file
    feature_counts_file = os.path.join(output_dir, f"feature_counts_{population}.tsv")
    feature_counts_df = pd.DataFrame(feature_counts, columns=['Column', 'Total Sum', 'Non-zero Count'])
    feature_counts_df.to_csv(feature_counts_file, sep='\t', index=False)

    # Save the dropped columns to a file
    dropped_columns_file = os.path.join(output_dir, f"dropped_columns_{population}.txt")
    with open(dropped_columns_file, 'w') as file:
        for col in dropped_columns:
            file.write(f"{col}\n")

    # Process the 'N2pair' column to keep NA as strings and integers as integers
    df['N2pair'] = df['N2pair'].apply(lambda x: 'NA' if pd.isna(x) else int(x) if str(x).isnumeric() else str(x))

    # Save the processed DataFrame to the output file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"\n\nProcessed file saved to: {output_file}")
    print(f"\nFeature counts saved to: {feature_counts_file}")
    print(f"\nDropped columns saved to: {dropped_columns_file}\n\n")

def main():
    parser = argparse.ArgumentParser(description="Process TSV file and save results.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input TSV file.")
    parser.add_argument("--output_file", type=str, required=True, help="Path where the processed TSV file will be saved.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory where the feature counts and dropped columns files will be saved.")
    parser.add_argument("--population", type=str, required=True, help="Population name to be used in output file names.")
    
    args = parser.parse_args()

    process_file(args.input_file, args.output_file, args.output_dir, args.population)

if __name__ == "__main__":
    main()

