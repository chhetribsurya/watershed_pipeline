import pandas as pd
import os
from os.path import join
import sys
import argparse

# Set pandas option to opt-in to future behavior
pd.set_option('future.no_silent_downcasting', True)

def process_data(out_dir, subjids_file, covariate_infile, pseudocount_infile, tpm_infile, population, tissue):
    """
    Processes data files for watershed model preps & analysis based on given parameters and generates output files.

    This function reads subject IDs, covariate files, pseudocount files, and TPM files,
    filters them based on the provided subject IDs, and writes the processed data to output files.

    Parameters:
    out_dir (str): The path to the output directory for output files.
    subjids_file (str): The path to the file containing subject IDs.
    covariate_infile (str): The path to the covariate input file.
    pseudocount_infile (str): The path to the pseudocount input file.
    tpm_infile (str): The path to the TPM input file.
    variable (str): The variable used to specify the data type (e.g., 'GLOBAL', 'EUR', 'AFR').
    tissue (str): The tissue type to be used in naming the output files.

    Example:
    process_data('/path/to/output_dir', 'subjids.csv', 'covariates.tab.gz', 'pseudocounts.tab', 'tpm.tab', 'GLOBAL', 'Liver')
    """
    # Ensure variable is one of the expected values
    #if variable not in ["GLOBAL", "EUR", "AFR"]:
    #    raise ValueError("Variable must be 'GLOBAL', 'EUR', or 'AFR'")

    # Define the peer_dir based on the variable
    #peer_dir = join(out_dir, "outputs", "1_expressions", f"PEER_{population}")
    #if not os.path.exists(peer_dir):
    #    os.makedirs(peer_dir)

    # Load SUBJID samples (reading only the first column)
    df_subjids = pd.read_csv(subjids_file, sep="\t", header=None, usecols=[0], names=['SUBJID'])
    print(f"\nTotal SUBJID samples loaded: {df_subjids.shape[0]}")

    # Check if the covariate file is gzipped and read accordingly
    if covariate_infile.endswith('.gz'):
        df_variable = pd.read_csv(covariate_infile, compression='gzip', sep="\t")
    else:
        df_variable = pd.read_csv(covariate_infile, sep="\t")

    # Process covariate file
    #df_variable = pd.read_csv(covariate_infile, compression="gzip", sep="\t")
    df_variable = df_variable[~df_variable.id.str.contains('PEER')]
    df_variable = df_variable.transpose()
    df_variable.columns = df_variable.iloc[0]
    df_variable = df_variable.drop('id')
    df_variable['sex'] = df_variable['sex'].replace({'XY': 1, 'XX': 2})
    cols = df_variable.columns.tolist()
    cols.remove('sex')
    cols.append('sex')
    df_variable = df_variable[cols]
    df_variable.reset_index(inplace=True)
    df_variable.rename(columns={'index': 'SUBJID'}, inplace=True)
    df_variable = df_variable[df_variable['SUBJID'].isin(df_subjids['SUBJID'])]
    print(f"Filtered {population} samples count: {df_variable.shape[0]}")
    #df_variable.to_csv(join(peer_dir, f'{tissue}_covariates.txt'), sep='\t', index=False)
    df_variable.to_csv(f'{tissue}_{population}_covariates.txt', sep='\t', index=False)

    # Process pseudocount file
    df_pseudocount = pd.read_csv(pseudocount_infile, sep="\t")
    df_pseudocount.rename(columns={df_pseudocount.columns[0]: 'Gene'}, inplace=True)
    for col in df_pseudocount.columns[1:]:
        df_pseudocount[col] = df_pseudocount[col].round().astype(int)
    df_pseudocount = df_pseudocount[['Gene'] + df_subjids['SUBJID'].tolist()]
    print(f"Filtered Read samples count: {len(df_subjids['SUBJID'].tolist())}")
    #df_pseudocount.to_csv(join(peer_dir, f'{tissue}_reads.txt'), sep='\t', index=False)
    df_pseudocount.to_csv(f'{tissue}_{population}_reads.txt', sep='\t', index=False)

    # Process TPM file
    df_tpm = pd.read_csv(tpm_infile, sep="\t")
    df_tpm.rename(columns={df_tpm.columns[0]: 'Gene'}, inplace=True)
    df_tpm = df_tpm[['Gene'] + df_subjids['SUBJID'].tolist()]
    print(f"Filtered TPM samples count: {len(df_subjids['SUBJID'].tolist())}")
    #df_tpm.to_csv(join(peer_dir, f'{tissue}_tpm.txt'), sep='\t', index=False)
    df_tpm.to_csv(f'{tissue}_{population}_tpm.txt', sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Process data files for watershed model preps & analysis.',
        usage='%(prog)s --out_dir OUT_DIR --subjids_file SUBJIDS_FILE --covariate_infile COVARIATE_INFILE '
              '--pseudocount_infile PSEUDOCOUNT_INFILE --tpm_infile TPM_INFILE --population VARIABLE --tissue TISSUE'
    )
    parser.add_argument('--out_dir', required=True, help='The path to the output directory for output files.')
    parser.add_argument('--subjids_file', required=True, help='The path to the file containing subject IDs.')
    parser.add_argument('--covariate_infile', required=True, help='The path to the covariate input file.')
    parser.add_argument('--pseudocount_infile', required=True, help='The path to the pseudocount input file.')
    parser.add_argument('--tpm_infile', required=True, help='The path to the TPM input file.')
    parser.add_argument('--population', required=True, help='The Population used to specify the data type (e.g., "GLOBAL", "EUR", "AFR").')
    parser.add_argument('--tissue', required=True, help='The tissue type to be used in naming the output files.')

    #if len(sys.argv) != 7:
    #    parser.print_help(sys.stderr)
    #    print("\n**Error: Incomplete set of arguments provided.\n")
    #    sys.exit(1)

    args = parser.parse_args()

    process_data(
        out_dir=args.out_dir,
        subjids_file=args.subjids_file,
        covariate_infile=args.covariate_infile,
        pseudocount_infile=args.pseudocount_infile,
        tpm_infile=args.tpm_infile,
        population=args.population,
        tissue=args.tissue
    )
