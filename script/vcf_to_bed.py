import pandas as pd
import numpy as np
import argparse
import gzip

def main(vcf, output):
    ###################
    # Read files
    ###################
    _MAX_RNAGE=5000
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

    def parse_info(info, keyword):
        if keyword in info: return info.split(keyword+'=')[1].split(';')[0]
        else: return np.nan

    for kw,dtype in zip(['AC','SVLEN','SVTYPE','STRAND','SUPPORT'],['float','float','str','str','float']):
        df_vcf[kw] = df_vcf["INFO"].apply(lambda x: parse_info(x, kw)).astype(dtype)
        
    df_vcf["_SVLEN"] = df_vcf["SVLEN"].fillna(1).abs()
    df_vcf.loc[df_vcf["_SVLEN"]>_MAX_RNAGE,'_SVLEN']=_MAX_RNAGE
    chromStart = df_vcf["POS"]-1
    chromEnd = df_vcf["POS"]+df_vcf["_SVLEN"]-1
    df_vcf["chromStart"] = chromStart
    df_vcf["chromEnd"] = chromEnd
    df_vcf["chromEnd"] = df_vcf["chromEnd"].astype(int)
    df_vcf[["CHROM","chromStart",'chromEnd']].to_csv(output, sep="\t", header=None, index=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="generate bed file from SV called vcf")
    parser.add_argument('--vcf')
    parser.add_argument('--output')
    args = parser.parse_args()
    
    main(args.vcf, args.output)