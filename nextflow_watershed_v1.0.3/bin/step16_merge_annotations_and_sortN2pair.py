#!/usr/bin/env python

import pandas as pd
import os
from os.path import join
import sys
import argparse
from sklearn.preprocessing import StandardScaler


def merge_annotation_files(annotation_dir, refbase_file, files_to_merge):
    """
    Merges multiple annotation files with a reference base file based on 'SubjectID' and 'GeneName'.

    Parameters:
    annotation_dir (str): Directory containing the annotation files.
    refbase_file (str): Path to the reference base file for merging.
    files_to_merge (list): List of annotation file names to be merged with the reference file.

    Returns:
    DataFrame: Merged DataFrame containing all annotations combined with the reference file.
    """

    # Read the reference base file
    refbase_df = pd.read_csv(refbase_file, sep='\t')
    print("Before drop reference base df:")
    print(refbase_df.head())

    # Drop the columns 'Zscore' and 'SignPVAL' - retain only PVAL column
    refbase_df.drop(columns=['Zscore', 'SignPVAL'], inplace=True)
    print("After drop reference base df:")
    print(refbase_df.head())

    # Drop the columns 'Zscore' and 'PVAL' - retain only SignPVAL column
    #refbase_df.drop(columns=['Zscore', 'PVAL'], inplace=True)
    #print("Reference base df after removal of cols:")
    #print(refbase_df.head())

    print(f"Original reference base file dimensions: {refbase_df.shape}")
    sys.stdout.flush()  
    print(refbase_df.head())
    sys.stdout.flush()  

    # Iteratively merge each file
    for file in files_to_merge:
        print(f"\n\nProcessing file: {file} ...\n")
        sys.stdout.flush()  
        
        file_df = pd.read_csv(join(annotation_dir, file), sep='\t')
        print(f"Original annotationfile dimensions: {file_df.shape}")
        sys.stdout.flush()  

        if not all(x in file_df.columns for x in ['SubjectID', 'GeneName']):
            raise ValueError(f"Missing required columns in {file}")

        # Merge using pandas merge
        refbase_df = pd.merge(refbase_df, file_df, on=['SubjectID', 'GeneName'], how='left')

        print(f"After merging {file} dimensions: {refbase_df.shape}\n")
        sys.stdout.flush()  
        print(refbase_df.head())
    
    return refbase_df


def merge_with_n2pairs_alt(merged_df, n2pair_df):
    """
    Merges the merged DataFrame with the N2 pairs DataFrame and rearranges the columns
    so that 'PVAL' is the second last and 'N2pair' is the last column.

    Parameters:
    merged_df (DataFrame): The DataFrame resulting from merging annotation files.
    n2pair_df (DataFrame): The DataFrame containing N2 pairs information.

    Returns:
    DataFrame: The final merged DataFrame with columns rearranged.
    """

    # Merge the merged_df with n2pair_df
    final_df = pd.merge(merged_df, n2pair_df, on=['SubjectID', 'GeneName'], how='left')

    # Rearrange columns so that PVAL is second last and N2pair is last
    if 'PVAL' in final_df.columns and 'N2pair' in final_df.columns:
        cols = list(final_df.columns)
        cols.remove('PVAL')
        cols.remove('N2pair')
        cols += ['PVAL', 'N2pair']
        final_df = final_df[cols]

    # if 'PVAL' in final_df.columns and 'N2pair' in final_df.columns:
       #  cols = [col for col in final_df.columns if col not in ['PVAL', 'N2pair']]
       #  cols += ['PVAL', 'N2pair']
       #  final_df = final_df[cols]

    return final_df



def merge_with_n2pairs(n2pair_df, merged_df):
    """
    Merges the merged DataFrame with the N2 pairs DataFrame and rearranges the columns
    so that 'PVAL' is the second last and 'N2pair' is the last column.

    Parameters:
    merged_df (DataFrame): The DataFrame resulting from merging annotation files.
    n2pair_df (DataFrame): The DataFrame containing N2 pairs information.

    Returns:
    DataFrame: The final merged DataFrame with columns rearranged.
    """

    # Merge the merged_df with n2pair_df
    final_df = pd.merge(n2pair_df, merged_df, on=['SubjectID', 'GeneName'], how='left')

    # Rearrange columns so that PVAL is second last and N2pair is last
    if 'PVAL' in final_df.columns and 'N2pair' in final_df.columns:
        cols = list(final_df.columns)
        cols.remove('PVAL')
        cols.remove('N2pair')
        cols += ['PVAL', 'N2pair']
        final_df = final_df[cols]

    return final_df


def sort_dataframe_by_n2pair(df):
    """
    Sorts the DataFrame such that NA values in 'N2pair' column are at the top and numeric values at the bottom.

    Parameters:
    df (DataFrame): The DataFrame to be sorted.

    Returns:
    DataFrame: Sorted DataFrame based on 'N2pair' column.
    """

    # Sorting such that NAs are at the top
    df_na = df[df['N2pair'].isna()]
    df_not_na = df[df['N2pair'].notna()]
    sorted_df = pd.concat([df_na, df_not_na])

    return sorted_df.reset_index(drop=True)


def main(args):
    # Set directories and annot files
    #pop = "GLOBAL" #args1
    #refbase_file = #args2 #join(datadir, "rv_expoutlier_baseref", pop, "gene-" + pop + "-rv.expoutlier.Pval.reference.tsv")
    #annotation_dir = #args3 "/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation"
    #n2pair_file = #args4 #"/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/GLOBAL/n2pair/gene-GLOBAL-rv.N2pairs.tsv"

    # Set parser arguments
    pop = args.pop
    annotation_dir = args.annotation_dir
    refbase_file = args.refbase_file
    n2pair_file = args.n2pair_file  

    #annotation_pop_dir = join(annotation_dir, pop, "annotation_inputfinal/genelevel_output")
    annotation_pop_dir = join(annotation_dir, pop)
    files_to_merge = [
          f"gene-{pop}-rv.CADD.collapse.tsv",
          f"gene-{pop}-rv.gencode.collapse.tsv",
          f"gene-{pop}-rv.afreq.collapse.tsv",
          f"gene-{pop}-rv.phyloP100way.collapse.tsv",
          f"gene-{pop}-rv.vep.loftee.collapse.tsv"
    ]

    # files_to_merge = [f"gene-{pop}-rv.{suffix}.tsv" for suffix in ["CADD", "gencode", "afreq", "phyloP100way", "vep.loftee"]]
    
    # Merging reference base with annotation files
    merged_df = merge_annotation_files(annotation_pop_dir, refbase_file, files_to_merge)

    # Read the n2pair file to combine with merged annotations # read as a nullable integer type
    n2pair_df = pd.read_csv(n2pair_file, sep="\t", dtype={"N2pair": pd.Int64Dtype()})

    # Merge dataframe with the N2 pairs
    final_merged_df = merge_with_n2pairs(n2pair_df, merged_df)

    # Sort such that NA values in 'N2pair' column are at the top
    final_sorted_df = sort_dataframe_by_n2pair(final_merged_df)
    
    # Replace NaN values in 'N2pair' column with 'NA'
    #final_sorted_df['N2pair'] = final_sorted_df['N2pair'].astype(str).replace('nan', 'NA')

    # Convert 'N2pair' column values to strings and replace 'nan' and '<NA>' with 'NA'
    final_sorted_df['N2pair'] = final_sorted_df['N2pair'].astype(str).replace({'nan': 'NA', '<NA>': 'NA'})
    
    # Save the final sorted DataFrame
    #final_output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.Zscorethresbased.tsv")
    final_output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.tsv")
    final_sorted_df.to_csv(final_output_file, sep='\t', index=False)
    print(f"\n Merged Annotation file with N2 pairs saved to: {final_output_file}")
    sys.stdout.flush()  
    
    # Save the annotations merged DataFrame
    #output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.tsv")
    #output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.Zscorethresbased.tsv")
    output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.Pvalthresbased.tsv")
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nMerged Annotation file saved to: {output_file}")
    sys.stdout.flush()  

    #########################################
    #                                       #
    # Standard scaling mean 0 std dev 1     #
    #                                       #   
    #########################################

    # Define columns to exclude from scaling
    exclude_columns = ['SubjectID', 'GeneName', 'N2pair'] + [col for col in final_sorted_df.columns if 'PVAL' in col]

    # Select columns to be scaled
    columns_to_scale = final_sorted_df.columns.difference(exclude_columns)

    # Apply standard scaling using sklearn
    scaler = StandardScaler()
    final_sorted_df[columns_to_scale] = scaler.fit_transform(final_sorted_df[columns_to_scale])
    #final_sorted_df[columns_to_scale] = scaler.fit_transform(final_sorted_df[columns_to_scale].fillna(0))

    # Save the final sorted DataFrame
    #final_output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.tsv")
    #final_output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.Zscorethresbased.stdscaled.tsv")
    final_output_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.tsv")
    final_sorted_df.to_csv(final_output_file, sep='\t', index=False)
    print(f"\n Merged Annotation file with N2 pairs saved to: {final_output_file}")
    sys.stdout.flush()  
    
    annotation_list = [
        "EncodeH3K36me3-max",  # Corrected hyphen
        "EncodeH3K79me2-max",  # Corrected hyphen
        "SIFTval",
        "downstream_gene_variant",
        "missense_variant",
        "LoF_LC",
        "SIFTcat_tolerated",
        "mamPhCons",
        "PolyPhenCat_benign",
        "upstream_gene_variant",
        "PolyPhenCat_possibly_damaging",
        "verPhCons",
        "PolyPhenCat_unknown",
        "SIFTcat_deleterious",
        "priPhCons",
        "PolyPhenVal",
        "EncodeH3K9me3-max",  # Corrected hyphen
        "PHRED",
        "splice_donor_variant",
        "EncodeH3K9ac-max",  # Corrected hyphen
        "EncodeH3K4me3-max",  # Corrected hyphen
        "EncodeH3K4me2-max"   # Corrected hyphen
    ]

    # Filter out the columns from final_sorted_df that are in annotation_list
    final_sorted_df_filtered = final_sorted_df.drop(columns=[col for col in annotation_list if col in final_sorted_df.columns])

    # Save the final sorted annotation filtered DataFrame
    final_output_filt_file = join(annotation_dir, pop, f"final-{pop}-rv.mergedannotation.plusN2pair.Pvalthresbased.stdscaled.annotfiltered.tsv")
    final_sorted_df_filtered.to_csv(final_output_filt_file, sep='\t', index=False)
    print(f"\n Merged filtered Annotation file with N2 pairs saved to: {final_output_file}")
    sys.stdout.flush()  
    
  
# if __name__ == "__main__":
#     main()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merges annotation files with a reference file and sorts the data based on N2pair column.")

    parser.add_argument("--pop", type=str, required=True, help="Population identifier used in file naming (e.g., GLOBAL, EUR, AFR).")
    parser.add_argument("--annotation_dir", type=str, required=True, help="Directory containing the annotation files.")
    parser.add_argument("--refbase_file", type=str, required=True, help="Path to the reference base file for merging.")
    parser.add_argument("--n2pair_file", type=str, required=True, help="Path to the N2 pairs file for merging.")

    args = parser.parse_args()
    main(args)

