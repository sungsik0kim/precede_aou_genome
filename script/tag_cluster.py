import pandas as pd
import numpy as np
import argparse
import gzip

def main(vcf, tag, output):
    ###################
    # Read files
    ###################
    skip_lines = 0
    with gzip.open(vcf, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_vcf = pd.read_csv(
        vcf, 
        sep='\t', 
        skiprows=skip_lines
    )
    df_vcf = df_vcf.rename(columns={'#CHROM': 'CHROM'})
    df_tag = pd.read_csv(tag, sep="\t", header=None)
    
    df_vcf["cluster_id"] = df_tag[3]
    df_vcf[['CHROM','POS','ID','cluster_id']].to_csv(output, index=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="tag cluster id to vcf")
    parser.add_argument('--vcf')
    parser.add_argument('--tag')
    parser.add_argument('--output')
    args = parser.parse_args()
    
    main(args.vcf, args.tag, args.output)