import os
import pandas as pd
import gzip
from tqdm import tqdm
import numpy as np
import argparse

def read_vcf_gz(vcf_gz):
    skip_lines = 0
    with gzip.open(vcf_gz, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_vcf = pd.read_csv(
        vcf_gz, 
        sep='\t', 
        skiprows=skip_lines
    )
    df_vcf = df_vcf.rename(columns={'#CHROM': 'CHROM'})
    return df_vcf

def main(vcf, vcf_unaffected, master, output):
    df_vcf = read_vcf_gz(vcf)
    df_unaffected = read_vcf_gz(vcf_unaffected)
    df_master = pd.read_csv(master)

    # Parse VCF table
    cols = df_master.loc[(~df_master['Analysis_pedigree_id'].isna()),'Analysis_sample_id'].tolist()
    cols = [col for col in cols if col in df_vcf.columns]
    df_vcf = df_vcf[["ID"]+cols].copy()
    m={'0/0':0 ,'0|0':0 ,'0/1':1 ,'0|1':1 ,'1|0':1 ,'1/1':1 ,'1|1':1, './.':'.'} # carrier status
    for col in cols:
        df_vcf[col] = df_vcf[col].str.split(':').str[0]
        df_vcf[col] = df_vcf[col].map(m)


    # Parse unaffected_rel_VCF table
    cols = df_master.loc[(~df_master['Analysis_pedigree_id'].isna()),'Analysis_sample_id'].tolist()
    cols = [col for col in cols if col in df_unaffected.columns]
    df_unaffected = df_unaffected[["ID"]+cols].copy()
    for col in cols:
        df_unaffected[col] = df_unaffected[col].str.split(':').str[0]
        df_unaffected[col] = df_unaffected[col].map(m)

    # Calculate Probability by Chance based on informative Meioses
    l_result=[]
    for _, row in df_vcf.iterrows():
        N=1 # Probability by Chance
        ID = row['ID']
        carriers = row.index[row==1]
        if any(carriers):
            mask_carrier = df_master['Analysis_sample_id'].apply(lambda x: x in carriers)
            pedigrees = df_master.loc[mask_carrier, 'Analysis_pedigree_id']
            pedigrees = pedigrees.unique().astype(int)

            for pedigree in pedigrees:
                mask_pedigree = df_master['Analysis_pedigree_id']==pedigree
                sample_id = df_master.loc[mask_pedigree, 'Analysis_sample_id'].tolist()
                sample_id=[sid for sid in sample_id if sid in df_vcf.columns.tolist()+df_unaffected.columns.tolist()]

                l_affected_status = [sid in df_vcf.columns.tolist() for sid in sample_id]
                l_carrier_status = []
                for sid, affected in zip(sample_id, l_affected_status):
                    if affected: carrier = df_vcf.loc[df_vcf["ID"]==ID, sid].iloc[0]
                    else: carrier = df_unaffected.loc[df_unaffected["ID"]==ID, sid].iloc[0]
                    l_carrier_status.append(carrier)
                m=0 # informative meosis
                for affected, carrier in zip(l_affected_status, l_carrier_status):
                    if affected==True and carrier==1: m+=1
                    elif affected==False and carrier==0: m+=1
                m = m-1 # Exclude index case
                _N = (1/2)**m
                N = N*_N
        l_result.append({
            'ID':ID
            ,'N':N
        })
    pd.DataFrame(l_result).to_csv(output, index=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="enrich pathogenic variants")
    parser.add_argument('--vcf')
    parser.add_argument('--unaffected_rel_vcf')
    parser.add_argument('--master')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(args.vcf, args.unaffected_rel_vcf, args.master, args.output)
