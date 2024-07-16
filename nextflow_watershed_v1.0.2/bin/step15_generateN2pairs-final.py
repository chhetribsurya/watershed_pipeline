#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from os.path import join
import sys


def process_rv_outlier_files(rv_file, exp_outlier_file):
    """
    Processes rare variant and expression outlier data.

    Parameters:
    rv_file (str): File path for the rare variant data.
    exp_outlier_file (str): File path for the expression outlier data.

    Returns:
    DataFrame: A DataFrame containing matched rare variants and expression outliers.
    """

    # Read the data files
    rv = pd.read_csv(rv_file, sep="\t")
    exp_outlier = pd.read_csv(exp_outlier_file, sep="\t")

    # Select genes that have at least one "outlier" status
    outlier_genes = exp_outlier.groupby('Gene').filter(lambda x: any(x['Status'] == 'outlier'))
    outlier_genes = outlier_genes.rename(columns={'Gene': 'GeneName', 'Individual': 'SubjectID'})

    # Count unique genes with outlier status
    outlier_uniq_genes_count = exp_outlier[exp_outlier['Status'] == 'outlier']['Gene'].nunique()
    print(f"\nTotal unique genes with at least one outlier individual status: {outlier_uniq_genes_count}\n")
    sys.stdout.flush()  

    # Rename columns for joining
    rv_rename = rv.rename(columns={'Gene': 'GeneName', 'Ind': 'SubjectID'})

    # Match rare variants with expression outliers
    outlier_rv_matched = pd.merge(rv_rename, outlier_genes, on=['GeneName', 'SubjectID'])

    # Reordering columns to make SubjectID the first column
    ordered_columns = ['SubjectID'] + [col for col in outlier_rv_matched.columns if col != 'SubjectID']
    outlier_rv_matched = outlier_rv_matched[ordered_columns]

    return outlier_rv_matched


def create_snp_alleles(df):
    """
    Adds a column to the DataFrame for SNP alleles.

    Parameters:
    df (DataFrame): The DataFrame to which the SNP alleles column will be added.

    Returns:
    DataFrame: The modified DataFrame with the SNP alleles column.
    """
    df['SNV_Alleles'] = df.Chrom + ':' + df.Start.astype(str) + ':' + df.Ref + ':' + df.Alt
    return df



def group_by_subjectid_genename(df):
    """
    Groups the DataFrame by SubjectID and GeneName and aggregates SNV alleles.

    Parameters:
    df (DataFrame): The DataFrame to be grouped.

    Returns:
    DataFrame: A DataFrame grouped by SubjectID and GeneName with aggregated SNV alleles.
    """
    return df.groupby(['SubjectID', 'GeneName']).agg(
        SNVs=('SNV_Alleles', lambda x: ','.join(np.sort(x.unique())))
    ).reset_index()



def group_by_subjectid_genename_alt(df):
    """
    Groups the DataFrame by 'SubjectID' and 'GeneName' and concatenates the 'SNV_Alleles' column.

    Parameters:
    df (pandas.DataFrame): The DataFrame containing 'SubjectID', 'GeneName', and 'SNV_Alleles' columns.

    Returns:
    pandas.DataFrame: A DataFrame with 'SubjectID', 'GeneName', and a new column 'SNVs', which
                      contains a comma-separated string of unique, sorted 'SNV_Alleles' for each group.
    """

    def concatenate_snvs(group):
        """Concatenate and sort unique SNV alleles into a comma-separated string."""
        sorted_snvs = sorted(group['SNV_Alleles'].unique())
        return ','.join(sorted_snvs)

    # Group by 'SubjectID' and 'GeneName', and apply the custom function to each group
    return df.groupby(['SubjectID', 'GeneName']).apply(concatenate_snvs).reset_index(name='SNVs')



def generate_n2_pairs(grouped_data):
    """
    Processes groups of the DataFrame, handling N2pair assignments and creating additional
    information for both N2 pairs and 'NA' cases.

    Parameters:
    grouped_data (DataFrame): The DataFrame with grouped data to be processed.

    Returns:
    DataFrame: Processed DataFrame with N2pair assignments.
    DataFrame: Additional DataFrame with extended information for N2 pairs, including SubjectIDs.
    DataFrame: Additional DataFrame with extended information for 'NA' cases, including SubjectIDs.
    """
    out = []
    extended_info_n2pairs = []
    extended_info_na = []
    n2_id = 0

    for (snvs, gene_name), group_df in grouped_data.groupby(['SNVs', 'GeneName']):
        group_df = group_df.copy()
        subject_ids = group_df['SubjectID'].tolist()  # Get list of SubjectIDs for the group

        if group_df.shape[0] > 1:
            # Handling N2 pairs
            n2pair = group_df.sample(n=2, random_state=32)
            n2pair['N2pair'] = n2_id
            group_df['N2pair'] = n2_id
            # Add information to extended_info_n2pairs
            extended_info_n2pairs.append({'N2pair': n2_id, 'SNVs': snvs, 'GeneName': gene_name,
                                          'RareVariantCount': snvs.count(',') + 1,
                                          'SubjectIDs': ';'.join(subject_ids)})
            n2_id += 1
            out.append(n2pair)
        else:
            # Handling 'NA' cases
            group_df['N2pair'] = 'NA'
            # Add information to extended_info_na
            extended_info_na.append({'N2pair': 'NA', 'SNVs': snvs, 'GeneName': gene_name,
                                     'RareVariantCount': snvs.count(',') + 1,
                                     'SubjectIDs': ';'.join(subject_ids)})
            out.append(group_df)

    # Creating DataFrames
    processed_df = pd.concat(out).reset_index(drop=True)
    extended_df_n2pairs = pd.DataFrame(extended_info_n2pairs)
    extended_df_na = pd.DataFrame(extended_info_na)

    return processed_df, extended_df_n2pairs, extended_df_na



def count_n2_pairs(df, output_file):
    """
    Analyzes the 'N2pair' column in the DataFrame to determine the count of N2 pairs and
    categorize elements into the 'watershed prediction dataset' and 'watershed training dataset'.
    The results are printed and saved to an output file.

    Parameters:
    df (DataFrame): The DataFrame containing the 'N2pair' column.
    output_file (str): Path to the output file where the results will be saved.

    Returns:
    None
    """

    # Count the total number of N2 pairs
    actual_pairings = df['N2pair'].nunique() - ('NA' in df['N2pair'].unique())

    # Number of elements in the 'watershed prediction dataset'
    prediction_dataset_count = df[df['N2pair'] != 'NA'].shape[0]

    # Number of elements in the 'watershed training dataset'
    training_dataset_count = df[df['N2pair'] == 'NA'].shape[0]

    # Prepare the results
    results = (
        f"Total number of N2 pairs (actual pairings): {actual_pairings}\n"
        f"Elements in the 'watershed prediction dataset': {prediction_dataset_count}\n"
        f"Elements in the 'watershed training dataset': {training_dataset_count}\n"
    )

    # Print and save the results
    print(results)
    with open(output_file, 'w') as file:
        file.write(results)



def main(args):

    # File paths
    # rv_file="/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/GLOBAL/gene-GLOBAL-rv.txt"
    # exp_outlier_file="/scratch16/abattle4/surya/datasets/WatershedAFR/data/data_prep/PEER_GLOBAL/kgpexV3.GLOBAL.outlier.controls.v3ciseQTLs_globalOutliersRemoved.txt"

    # File paths
    rare_variant_file = args.rare_variant_file
    exp_outlier_file = args.exp_outlier_file
    outfile = args.outfile

    # Load and process rare variant and expression outlier data
    outlier_rv_matched = process_rv_outlier_files(rare_variant_file, exp_outlier_file)
    snv_annotation = create_snp_alleles(outlier_rv_matched)
    print("\n SNV variant IDs generated ...\n")
    sys.stdout.flush()  

    # Grouping data by subjectID and gene
    grouped_data = group_by_subjectid_genename(snv_annotation)
    grouped_data1 = group_by_subjectid_genename_alt(snv_annotation)
    print("\n Grouping by Indiv-Gene completed ...\n")
    sys.stdout.flush()  

    # Asserting equality of the dataframes based on alternative approach
    assert grouped_data1.equals(grouped_data), "The dataframes are not equal"
    print("\n The dataframes are equal based on alt. method \n")
    sys.stdout.flush()  

    # Clustering/grouping based on Gene and SNV IDs, i.e., matching all rvs for the same gene
    processed_df, extended_df_n2pairs, extended_df_na = generate_n2_pairs(grouped_data)
    print("\n Grouping/Clustering by gene-and-SNVvariantID completed ...\n")
    sys.stdout.flush()  
   
    # Annotate filenames 
    processed_df_file = outfile.replace('.N2pairs.tsv', '.n2.na.info.txt')
    extended_df_n2pairs_file = outfile.replace('.N2pairs.tsv', '.n2.info.txt')
    extended_df_na_file = outfile.replace('.N2pairs.tsv', '.na.info.txt')

    # Save above files
    processed_df.to_csv(processed_df_file, sep='\t', index=False, header=True) 
    extended_df_n2pairs.to_csv(extended_df_n2pairs_file, sep='\t', index=False, header=True)
    extended_df_na.to_csv(extended_df_na_file, sep='\t', index=False, header=True)

    # Final N2pair dataframe
    train_na_test_n2pair_df = processed_df.drop('SNVs', axis=1)
    train_na_test_n2pair_df.to_csv(outfile, sep="\t", index=False, header=True)
    print(f"\n N2pair annotation file saved to: {outfile} \n")
    sys.stdout.flush()  

    # Counting N2 pairs and saving to file
    n2_pairs_file = outfile.replace('.N2pairs.tsv', '.n2.na.stats.txt')
    count_n2_pairs(train_na_test_n2pair_df, n2_pairs_file)
    print(f"\n N2pairs statistics saved to: {n2_pairs_file} \n") 
    sys.stdout.flush()  
    print("\n N2pairs generation completed ...\n")
    sys.stdout.flush()  


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and analyze N2pair annotation data.")
    parser.add_argument("--rare_variant_file", type=str, required=True, help="Rare variants file path")
    parser.add_argument("--exp_outlier_file", type=str, required=True, help="Expression outlier file path")
    parser.add_argument("--outfile", type=str, required=True, help="File path for the final result and N2 pair annotations")
    
    args = parser.parse_args()
    main(args)


#python generateN2pairs-final.py --rare_variant_file path/to/rare_variant_file.tsv --exp_outlier_file path/to/exp_outlier_file.tsv --outfile path/to/output_file.tsv

