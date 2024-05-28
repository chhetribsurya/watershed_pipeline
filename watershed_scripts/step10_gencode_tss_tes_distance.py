#!/usr/bin/env python

# conda activate py3-env

"""
Converts given IPYNB to a Python Script with Comments, Details, and Usage Help.
Script Name: script_name.py
Usage: python script_name.py

This script obtains the distance from a rare variant to the TSS or TES. It calculates:
- distTSS: Absolute distance to the transcription start site.
- distTES: Absolute distance to the transcription end site.
"""

import os
import numpy as np
import pandas as pd

def infostr_to_dict(infostr, kv_sep=' ', sep=';'):
    """
    Converts a string to a dictionary.

    Parameters:
    ----------
    infostr : str
        Info column from gencode.
    
    kv_sep : str
        Separator between key and value.
    
    sep : str
        Separator between key-value pairs.

    Returns:
    -------
    dict
        Dictionary representation of the input string.
    """
    
    kv_pairs = infostr.rstrip(sep).split(sep)

    infodict = {}
    for pair in kv_pairs:
        pair = pair.lstrip()

        k, v = pair.split(kv_sep)
        v = v.replace('"','')
        infodict[k] = v
    
    return infodict

def infostr_to_gene_id(infostr):
    """
    Extracts gene_id from column 9 in gencode.

    Parameters:
    ----------
    infostr : str
        Info column from gencode.

    Returns:
    -------
    str
        Extracted gene_id.
    """
    return infostr_to_dict(infostr)['gene_id']

def get_distance(rv_row):
    """
    Returns a rare variant's absolute distance to TSS and TES.

    Parameters:
    ----------
    rv_row : pandas.Series
        Row from rare variant file as DataFrame.

    Returns:
    -------
    tuple
        distTSS : int
            Absolute distance to transcription start site.
        distTES : int
            Absolute distance to transcription end site.
    """
    
    pos = rv_row['Start'] + 1
    strand = rv_row['Strand']
    
    if strand == '+':
        distTSS = np.abs(rv_row['genStart'] - pos)
        distTES = np.abs(rv_row['genEnd'] - pos)
    else:
        distTSS = np.abs(rv_row['genEnd'] - pos)
        distTES = np.abs(rv_row['genStart'] - pos)
    
    return distTSS, distTES

def main():

    #gencode_file = '/scratch16/abattle4/victor/victor/WatershedAFR/raw_data/GTEx/gencode.v26.GRCh38.genes.gtf'
    gencode_file="/scratch16/abattle4/surya/datasets/WatershedAFR/raw_data/1KG/gencode.v38.GRCh38.genes.gtf"
    
    #rv_file = '/scratch16/abattle4/victor/victor/WatershedAFR/data/rare_variants_gnomad/gene-AFR-rv.txt'
    rv_file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt'
    
    #gencode = pd.read_csv(gencode_file, names=['Chrom', 'Source', 'Feature', 'genStart', 'genEnd', 'Score','Strand', 'Phase', 'Info'], sep='\t', index_col=False, skiprows=6)
    colnames = ['Chrom', 'Source', 'Feature', 'genStart', 'genEnd', 'Score','Strand', 'Phase', 'Info']
    gencode = pd.read_csv(gencode_file, names=colnames, sep='\t', index_col=False, skiprows=5)
    
    #rv = pd.read_csv(rv_file, sep='\t', index_col=False)
    rv = pd.read_csv(rv_file, sep='\t', index_col=False)
    
    # Parse gene ID and merge with rare variant file
    gencode['Gene'] = gencode['Info'].apply(infostr_to_gene_id)
    gencode_gene = gencode[gencode['Feature'] == 'gene'][['Gene', 'genStart', 'genEnd', 'Strand']]
    gencode_anno = pd.merge(rv, gencode_gene, how='left', on='Gene')
    
    # GET TSS TES distances
    gencode_anno['distTSS'], gencode_anno['distTES'] = zip(*gencode_anno.apply(get_distance, axis=1))
    gencode_anno_df = gencode_anno[['Gene', 'Ind', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF', 'distTSS', 'distTES']]
    
    # Output file
    gencode_anno_file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/gencode/gene-GLOBAL-rv.gencode.txt'
    gencode_anno_df.to_csv(gencode_anno_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
