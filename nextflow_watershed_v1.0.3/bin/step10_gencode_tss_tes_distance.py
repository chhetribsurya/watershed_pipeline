#!/usr/bin/env python

# conda activate py3-env

"""
This script obtains the distance from a rare variant to the TSS or TES. It calculates:
- distTSS: Absolute distance to the transcription start site.
- distTES: Absolute distance to the transcription end site.
"""

import os
import numpy as np
import pandas as pd
import sys

#gencode_infile = sys.argv[1]
#rv_infile = sys.argv[2]
#outfile = sys.argv[3]

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

def main(gencode_file, rv_file, output_file):
 
    colnames = ['Chrom', 'Source', 'Feature', 'genStart', 'genEnd', 'Score','Strand', 'Phase', 'Info']
    
    # Load gencode file
    gencode = pd.read_csv(gencode_file, names=colnames, sep='\t', index_col=False, skiprows=5)
   
    # Load rarevariant file 
    rv = pd.read_csv(rv_file, sep='\t', index_col=False)
    
    # Parse gene ID and merge with rare variant file
    gencode['Gene'] = gencode['Info'].apply(infostr_to_gene_id)
    gencode_gene = gencode[gencode['Feature'] == 'gene'][['Gene', 'genStart', 'genEnd', 'Strand']]
    gencode_anno = pd.merge(rv, gencode_gene, how='left', on='Gene')
    
    # GET TSS TES distances
    gencode_anno['distTSS'], gencode_anno['distTES'] = zip(*gencode_anno.apply(get_distance, axis=1))
    gencode_anno_df = gencode_anno[['Gene', 'Ind', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'AF', 'distTSS', 'distTES']]
    
    # Output file
    gencode_anno_file = output_file
    gencode_anno_df.to_csv(gencode_anno_file, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python step10_gencode_tss_tes_distance.py <GENCODE_FILE> <RV_FILE> <OUTPUT_FILE>")
        sys.exit(1)

    gencode_file = sys.argv[1]
    rv_file = sys.argv[2]
    output_file = sys.argv[3]
    main(gencode_file, rv_file, output_file)


