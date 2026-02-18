import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy.stats import pearsonr
import argparse

def variant_merge(vep, af, filt, cluster, cadd_sv, output):
    ###################
    # Read files
    ###################
    df_vep = pd.read_csv(vep)
    print(df_vep.value_counts("SVTYPE"))
    skip_lines = 0
    with open(af, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_af = pd.read_csv(af)
    df_mask = pd.read_csv(filt)
    df_cluster = pd.read_csv(cluster)
    df_caddsv = pd.read_csv(cadd_sv, sep="\t")
    sample_cols = df_mask.columns[3:].tolist()
    
    ## Masking low confident loci and calls
    df_merge = pd.merge(df_mask[['ID']], df_vep, on="ID",how="left")
    for sample in tqdm(sample_cols):
        gt = df_merge[sample].values
        gt[~df_mask[sample]]=np.nan
        df_merge[sample] = gt
    print(df_merge.value_counts("SVTYPE"))
    
    ## Join AF
    af_cols = df_af.columns.tolist()[-14:]
    df_merge = pd.merge(df_merge, df_af[['ID']+af_cols], on="ID",how="left")
    print(df_merge.value_counts("SVTYPE"))
    ## Join cluster tag
    df_merge = pd.merge(df_merge, df_cluster, on=['CHROM','POS','ID'],how="left")
    ## Join CADD-SV
    df_caddsv = df_caddsv[['Name','CADD-SV_PHRED-score']]
    df_caddsv.columns = ['ID','CADD-SV_PHRED-score']
    df_merge = pd.merge(df_merge, df_caddsv, on='ID',how="left")
    
    ## Save
    df_merge.to_csv(output, index=False)
    print(df_merge.value_counts("SVTYPE"))
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Merge and Filter Variants")
    parser.add_argument('--vep')
    parser.add_argument('--af')
    parser.add_argument('--filter')
    parser.add_argument('--cluster')
    parser.add_argument('--cadd_sv')
    parser.add_argument('--output')
    args = parser.parse_args()
    
    variant_merge(args.vep, args.af, args.filter, args.cluster, args.cadd_sv, args.output)