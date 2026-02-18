import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import argparse


def vep_parse(vep_path, output):
    ###################
    # ## Read vep data
    ###################
    skip_lines = 0
    with open(vep_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df = pd.read_csv(
        vep_path, 
        sep='\t', 
        skiprows=skip_lines
    )
    df = df.rename(columns={'#CHROM': 'CHROM'})

    sample_cols= df.columns[9:].tolist()
    aou_cols = [col for col in sample_cols if 'sniffles2.hg38' in col]
    cohort_cols = list(set(sample_cols)-set(aou_cols))

    ###################
    # ## Parse df
    ###################
    df = df.rename(columns={'#CHROM': 'CHROM'})

    def parse_info(info, keyword):
        if keyword in info: return info.split(keyword+'=')[1].split(';')[0]
        else: return np.nan

    for kw,dtype in zip(['AC','SVLEN','SVTYPE','STRAND','SUPPORT'],['float','float','str','str','float']):
        df[kw] = df["INFO"].apply(lambda x: parse_info(x, kw)).astype(dtype)

    vep_cols = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|GIVEN_REF|USED_REF|BAM_EDIT|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|SV_overlap_AF|SV_overlap_PC|SV_overlap_name|PHENOTYPES|pLI_gene_value|OpenTargets_geneId|OpenTargets_l2g'.split('|')
    for i, kw in tqdm(enumerate(vep_cols)):
        if kw=='STRAND': continue
        df[kw] = df['INFO'].str.split('CSQ=').str[1].str.split('|').str[i]

    
    ###################
    # ## Filtering
    ###################
    #### Variant filtering step 2: remove low GQ & high missing GT rate
    _df_gt={}
    for col in tqdm(sample_cols):
        _df_gt[col]=df[col].str.split(':').str[0]
    df_gt = pd.DataFrame(_df_gt)[sample_cols]

    d = {'./.':np.nan, '0/0':0, '0/1':1, '1/1':2}
    for col in tqdm(sample_cols):
        df_gt[col] = df_gt[col].map(d)

    for col in tqdm(sample_cols):
        df[col] = df_gt[col]
    df.to_csv(output, index=False)
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Filter and Parse VEP result")
    parser.add_argument('--vep', help="vep result in vcf format")
    parser.add_argument('--output', help="output file")
    arags = parser.parse_args()
    
    vep_parse(arags.vep, arags.output)
    
    
    
    
    