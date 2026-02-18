import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import os


def main(plink_ld, outcsv):
    df = pd.read_csv(plink_ld, sep="\t").rename(columns={'#CHROM_A':'CHROM_A'})
    df = df.loc[(df["ID_B"].str.split(':').str[2].str.len()<=2) & (df["ID_B"].str.split(':').str[3].str.len()<=2)] # limit indel size
    df["dist"] = (df['POS_B']-df['POS_A']).abs()
    df = df.loc[df['dist']>50] # remove potential snp inside sv
    
    df_sorted = df.sort_values(by=['ID_A', 'UNPHASED_R2','dist'], ascending=[True, False, True])
    df_subset = df_sorted.groupby('ID_A').head(3) # get max 3 snps per sv
    df_subset = df_subset.reset_index(drop=True)
    df_subset['ID_B_2'] = df_subset['ID_B'].apply(lambda x:'-'.join(x.split(':')).replace('chr',''))
    df_subset.to_csv(outcsv, sep='\t', index=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--plink_ld')
    parser.add_argument('--outcsv')
    args = parser.parse_args()
    main(args.plink_ld, args.outcsv)

