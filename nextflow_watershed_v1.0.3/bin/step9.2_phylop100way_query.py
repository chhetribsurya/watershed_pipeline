import os
import argparse

import pandas as pd
import pysam


def query_phylop(chrom, pos):
    """Queries the phylop score at the rare variant

    Parameters
    ----------
    chrom : str
      Chromosome the rare variant is located on
    pos : int
      0-based position of the rare variant

    Returns
    -------
    str
    """
    try:
        query_score = next(phylop.fetch(chrom, pos, pos + 1, parser=pysam.asTuple()))[3]
    except StopIteration:
        query_score = "0"  # impute as 0
    return query_score

parser = argparse.ArgumentParser()
parser.add_argument(
    "--phylop", required=True, help="bed file containing the phylop 100 way scores"
)
parser.add_argument(
    "--rarevariant", required=True, help="bed file containing the rare variants"
)
parser.add_argument(
    "--pop", required=True, help="string containing name of population"
)
parser.add_argument(
    "--outdir", required=True, help="string containing the path for output directory"
)
args = parser.parse_args()

phylop_file = args.phylop
rv_file = args.rarevariant
pop_name = args.pop
outputdir = args.outdir

# create output dir
os.makedirs(outputdir, exist_ok=True)

# output file location
out_file = os.path.join(outputdir, "gene-" + pop_name + "-rv.phyloP100way.bed")

# open PhyloP 100 way scores from bed file
phylop = pysam.TabixFile(phylop_file)

out_writer = open(out_file, "w", buffering=1)

with open(rv_file, "r") as f:
    for idx, line in enumerate(f):
        chrom, start, end = line.rstrip().split("\t")
        score = query_phylop(chrom, int(start))
        out_writer.write("\t".join([chrom, start, end, score]) + "\n")

out_writer.close()
print(out_file)
