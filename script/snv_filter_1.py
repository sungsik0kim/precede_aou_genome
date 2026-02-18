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

def main(vcf, case, ctrl, outdir):
    ###################
    # Read VCF
    ###################
    filt_variant_bed = os.path.join(outdir,'filt_variant.bed')
    filt_variant_id = os.path.join(outdir,'filt_variant.id')
    filt_sample =  os.path.join(outdir,'filt_sample.txt')
    
    df_vcf = read_vcf_gz(vcf)
    df_ctrl =pd.read_csv(ctrl)
    df_case =pd.read_csv(case)
    
    sample_case = [sample for sample in  df_case['analysis_sample_id'].astype(str) if sample in df_vcf.columns.tolist()]
    sample_ctrl = [sample for sample in  df_ctrl['analysis_sample_id'].astype(str) if sample in df_vcf.columns.tolist()]

    
    df_vcf_case = df_vcf[sample_case].copy()
    mask_non_ref = []
    for i in tqdm(range(df_vcf_case.shape[0])):
        uniq = df_vcf_case.loc[i].str.split(':').str[0].unique()
        uniq = set(uniq)
        if len(uniq - set(['./.', '0/0']))>0: mask_non_ref.append(True)
        else: mask_non_ref.append(False)
    mask_PASS = (df_vcf['FILTER']=='PASS').tolist()
    mask_non_ref = np.array(mask_non_ref)
    mask_PASS = np.array(mask_PASS)

    print('Num Original Variants: {}'.format(df_vcf.shape[0]))
    print('Variants in Case: {}'.format(df_vcf.loc[mask_non_ref].shape[0]))
    print('Num FILTER=PASS Variants in Case: {}'.format(df_vcf.loc[mask_non_ref & mask_PASS].shape[0]))
    
    df_vcf_filt = df_vcf.loc[mask_non_ref & mask_PASS].copy()
    df_vcf_filt_bed = df_vcf_filt[['CHROM','POS','POS']].copy()
    df_vcf_filt_bed.columns = ['CHROM','START','END']
    df_vcf_filt_bed['START'] = df_vcf_filt_bed['START']-1
    df_vcf_filt_bed.to_csv(filt_variant_bed, sep="\t", index=False, header=None)

    with open(filt_sample, 'w') as f:
        for s in sample_case: f.write(s+'\n')
        for s in sample_ctrl: f.write(s+'\n')
    with open(filt_variant_id, 'w') as f:
        for ID in df_vcf_filt["ID"]: f.write(ID+'\n')
    
            
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Filter VCF")
    parser.add_argument('--vcf', help="vcf to be filtered")
    parser.add_argument('--case', help="csv for case")
    parser.add_argument('--ctrl', help="csv for ctrl")
    parser.add_argument('--outdir')
    args = parser.parse_args()
    
    main(args.vcf, args.case, args.ctrl, args.outdir)
