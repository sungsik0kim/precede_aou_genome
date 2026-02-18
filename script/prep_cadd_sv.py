import argparse
import gzip
import pandas as pd

def expand_col(series, tag, dtype=None):
    if dtype:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0].astype(dtype)
    else:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0]


def main(vcf_in,bed_out):
    print('###### Read and Parse Data #######')
    skip_lines = 0
    with gzip.open(vcf_in, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df = pd.read_csv(
        vcf_in, 
        sep='\t', 
        skiprows=skip_lines
    )
    df = df.rename(columns={'#CHROM': 'CHROM'})
    
    df['SVTYPE'] = expand_col(df['INFO'], 'SVTYPE')
    df['SVLEN'] = expand_col(df['INFO'], 'SVLEN').astype(float)
    
    df_bed = df.loc[df['SVTYPE']!="BND"].copy()
    df_bed["END"] = df["POS"]+df["SVLEN"].abs()
    df_bed["END"] = df_bed["END"].astype(int)
    df_bed = df_bed[["CHROM",'POS','END',"SVTYPE","ID"]]
    df_bed["POS"]=df_bed["POS"]-1
    df_bed.to_csv(bed_out, header=None, index=False, sep="\t")


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Prep CADD-SV input")
    parser.add_argument('--vcf_in')
    parser.add_argument('--bed_out')
    args = parser.parse_args()
    main(args.vcf_in, args.bed_out)

