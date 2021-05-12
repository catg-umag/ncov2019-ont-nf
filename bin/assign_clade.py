#!/usr/bin/env python3

import sys
import io
import os
import pandas as pd
import gzip

def read_vcf(path):
    with gzip.open(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def main(argv):
    df = read_vcf(argv[0])
    clades = pd.read_csv(argv[1],sep="\t")
    assign_clade = ''
    for clade in ["S","L","V","GH","GR","GV","GRY","G"]:
       variants = clades.loc[clades.CLADE == clade]
       idx_cols = ["POS","REF","ALT"] if clade != "L" else ["POS","REF"]
       isin = variants.set_index(idx_cols).index.isin(df.set_index(idx_cols).index).all()    # if clade L, if clade GRY
       if(isin):
          assign_clade = clade
          break
    print(assign_clade)
    return assign_clade

if __name__ == "__main__":
    main(sys.argv[1:])
