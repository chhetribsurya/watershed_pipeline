#!/usr/bin/env python

import os
import argparse
import glob
import pandas as pd
import pysam


def parse_vcf(vcf_path, anno_list):
    vcf = pysam.VariantFile(vcf_path)

    records = []
    for rec in vcf.fetch():
        info_data = rec.info['CSQ']
        # Initialize dictionary to hold annotation flags
        annotation_flags = {anno: 0 for anno in anno_list}
        
        for info in info_data:
            fields = info.split('|')
            consequences = fields[1].split('&')
            lof = fields[-4]  # Assuming LoF is the fourth last field in CSQ
            
            # Update flags for consequences
            for consequence in consequences:
                if consequence in annotation_flags:
                    annotation_flags[consequence] = 1

            # Update flags for LoF
            if lof == 'HC':
                annotation_flags['LoF_HC'] = 1
            elif lof == 'LC':
                annotation_flags['LoF_LC'] = 1

        record = {
            'Chrom': rec.chrom,
            'Pos': rec.pos,
            'Ref': rec.ref,
            'Alt': ','.join([alt for alt in rec.alts]),
            **annotation_flags
        }
        records.append(record)

    return pd.DataFrame(records)


def process_annotations(df):
    # Collapse transcripts by taking maximum to get SNV level annotations
    snv_df = df.groupby(["Chrom", "Pos", "Ref", "Alt"], as_index=False).max()

    # Set LoF_LC to 0 if LoF_HC is 1
    snv_df["LoF_LC"] = [0 if hc == 1 else lc for hc, lc in zip(snv_df["LoF_HC"], snv_df["LoF_LC"])]

    return snv_df


def main():
    usage = "Parses VCF output from VEP/LOFTEE to a table where rows are rare variants and columns are the annotations"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument(
        "--anno_vcf_list",
        required=True,
        help="VCF file list  of annotations from VEP and LOFTEE in vcf.gz format.",
    )
    parser.add_argument(
        "--outputfile",
        required=True,
        help="Output file path for the combined SNV table.",
    )

    args = parser.parse_args()
    anno_vcf_files = args.anno_vcf_list
    outputfile = args.outputfile
 
    # List of annotations to extract
    anno_list = [
        "3_prime_UTR_variant",
        "5_prime_UTR_variant",
        "TF_binding_site_variant",
        "downstream_gene_variant",
        "intergenic_variant",
        "intron_variant",
        "missense_variant",
        "non_coding_transcript_exon_variant",
        "regulatory_region_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
        "stop_gained",
        "synonymous_variant",
        "upstream_gene_variant",
        "LoF_HC",
        "LoF_LC",
    ]

    # Process and combine dataframes for each chromosome
    anno_snv_all_df = None
    #for vcf_file in list(anno_vcf_files):
    for vcf_file in anno_vcf_files.split(","):
        #if not os.path.exists(vcf_file):
        #    print(f"File not found: {vcf_file}")
        #    continue

        # Print the VCF file being processed
        print(f"Processing VEP file parse: {vcf_file}")

        # Parse the VCF and process the annotations
        df = parse_vcf(vcf_file, anno_list)
        anno_snv_df = process_annotations(df)

        if anno_snv_all_df is None:
            anno_snv_all_df = anno_snv_df
        else:
            anno_snv_all_df = pd.concat([anno_snv_all_df, anno_snv_df])

    # Save the combined dataframe
    anno_snv_all_df.to_csv(outputfile, sep="\t", index=False)
    print(f"\n**** SNV level VEP and LOFTEE annotations are saved to:\n{outputfile}\n")

if __name__ == "__main__":
    main()

