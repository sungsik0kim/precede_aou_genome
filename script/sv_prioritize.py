import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import os

def expand_col(series, tag, dtype=None):
    if dtype:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0].astype(dtype)
    else:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0]
    
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

def read_vcf(vcf):
    skip_lines = 0
    with open(vcf, 'rt') as f:
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
    return df_vcf

def main(filt2_csv, af, vep, cadd, annotSV, family_segregation, master, outcsv):
    ### Read Data ###
    df_vep = read_vcf(vep)
    df_cadd = pd.read_csv(cadd, sep='\t').rename(columns={'#Chrom': 'Chrom'})
    df_annotsv = pd.read_csv(annotSV, sep='\t')
    df_segregation = pd.read_csv(family_segregation)
    df_filt2 = pd.read_csv(filt2_csv)
    df_af = pd.read_csv(af)
    df_segregation = pd.read_csv(family_segregation)
    df_master = pd.read_csv(master)
    
    #### Parse VEP ###
    df_vep["_vep"] = expand_col(df_vep['INFO'], 'CSQ')
    vep_cols='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|GIVEN_REF|USED_REF|BAM_EDIT|SOURCE|HGVS_OFFSET|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|ClinVar|ClinVar_CLNSIG|ClinVar_CLNDN|phyloP100way|phyloP100way_min|phyloP100way_mean|phyloP100way_max|phastCons100way|phastCons100way_min|phastCons100way_mean|phastCons100way_max|GERP|GERP_min|GERP_mean|GERP_max'.split('|')
    for enum, col in enumerate(vep_cols):
        df_vep[col]= df_vep["_vep"].str.split('|').str[enum]
    vep_slim_cols = [
                        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                        'HGVSc','HGVSp','ClinVar', 'ClinVar_CLNSIG', 'ClinVar_CLNDN',
                        'phyloP100way_min', 'phyloP100way_mean', 'phyloP100way_max',
                        'phastCons100way_min','phastCons100way_mean', 'phastCons100way_max',
                        'GERP_min', 'GERP_mean', 'GERP_max'                 
                    ] 
    df_vep_slim = df_vep[vep_slim_cols].copy()
    cols = [
            'phyloP100way_min', 'phyloP100way_mean', 'phyloP100way_max',
            'phastCons100way_min','phastCons100way_mean', 'phastCons100way_max',
            'GERP_min', 'GERP_mean', 'GERP_max'                 
           ]
    for col in cols:
        df_vep_slim.loc[df_vep_slim[col]=='',col] = np.nan
        df_vep_slim[col] = df_vep_slim[col].astype(float)
    
    
    #### Parse AnnotSV ###
    df_annotsv_full = df_annotsv.loc[df_annotsv['Annotation_mode']=='full']
    df_annotsv_split = df_annotsv.loc[df_annotsv['Annotation_mode']=='split']
    annotsv_cols = ['AnnotSV_ID',
                    'SV_length',
                    'SV_type',
                    'Annotation_mode',
                    'Gene_name',
                    'Closest_left',
                    'Closest_right',
                    'Tx',
                    'RE_gene',
                    'ACMG',
                    'HI',
                    'TS',
                    'GnomAD_pLI',
                    'ExAC_pLI',
                    'AnnotSV_ranking_score',
                    'AnnotSV_ranking_criteria',
                    'ACMG_class']
    df_annotsv_full_slim = df_annotsv_full[annotsv_cols]

    # -- Add info from "Split" data frame
    l_res=[] 
    for ID, df_gb in df_annotsv_split[["AnnotSV_ID","Tx","Location","Dist_nearest_SS"]].groupby('AnnotSV_ID'):
        l_res.append({
            'AnnotSV_ID': ID
            ,'Tx': ';'.join(df_gb["Tx"])
            ,'Location': ';'.join(df_gb["Location"])
            ,'Dist_nearest_SS': ';'.join(df_gb["Dist_nearest_SS"].astype(str))

        })
    df_annotsv_split2 = pd.DataFrame(l_res)

    df_annotsv_full_slim = pd.merge(df_annotsv_full_slim, df_annotsv_split2, on='AnnotSV_ID', how="left")


    df_annotsv_full_slim.columns = ["AnnotSV_"+col for col in df_annotsv_full_slim.columns]
    df_annotsv_full_slim = pd.merge(df_annotsv_full_slim, df_annotsv_full[["AnnotSV_ID",'ID']]
                                    , left_on="AnnotSV_AnnotSV_ID", right_on="AnnotSV_ID", how='left')
    df_annotsv_full_slim = df_annotsv_full_slim.drop(columns=['AnnotSV_AnnotSV_ID'])
    
    ### CADD-SV ###
    cadd_slim_cols = [
        'CADD-SV_PHRED-score',
        'CADD-SV_Raw-score',
    ]
    df_cadd_slim = df_cadd[cadd_slim_cols]
    df_cadd_slim["ID"] = df_cadd["Name"].copy()
    
    ### Family Segregation ###
    df_segregation["CHROM"]=df_segregation["ID"].str.split('_').str[0]
    df_segregation["POS"]=df_segregation["ID"].str.split('_').str[1]
    df_segregation["REF"]=df_segregation["ID"].str.split('_').str[2]
    df_segregation["ALT"]=df_segregation["ID"].str.split('_').str[3]
    df_segregation["Prob_random_segregation"]=df_segregation["N"]
    
    
    ### Merge ###
    df_merge = pd.merge(df_vep_slim, df_annotsv_full_slim, on="ID", how="left")
    df_merge = pd.merge(df_merge, df_cadd_slim, on="ID", how="left")
    df_merge = pd.merge(df_merge, df_af[['ID','gnomAD_Max_PopMax_AF','CoLoRdb_AF','HPRC_AF']], on="ID", how="left")
    df_merge = pd.merge(df_merge, df_filt2[['ID','case_EAF', 'ctrl_AF','AoU_Max_AF','internal_ctrl_AF'
                                            ,'Unaffected_Rel_EAF', 'OR', 'FPC_sample_id']], on="ID", how="left")
    df_merge = pd.merge(df_merge, df_segregation[['ID','Prob_random_segregation']], on="ID", how="left")
    
    l_ambry_accession=[]
    for sample_ids in df_merge["FPC_sample_id"]:
        sample_ids = sample_ids.split(', ')
        _l_ambry_accession=[]
        for sample_id in sample_ids:
            _l_ambry_accession.append(str(df_master.loc[df_master['Analysis_sample_id']==sample_id,'Ambry accession #'].iloc[0]))
        l_ambry_accession.append(', '.join(_l_ambry_accession))
    df_merge['Ambry_accession'] = l_ambry_accession
    df_merge['Ambry_accession'] = df_merge['Ambry_accession'].str.replace('nan','N/A')
    
    ### Prioritize & Print ###
    mask1 = df_merge['ClinVar_CLNSIG'].astype(str).apply(lambda x: ("pathogenic" in x) or ("Pathogenic" in x) )
    mask2 = df_merge['phyloP100way_max']>5
    mask3 = (df_merge['phastCons100way_max']>0.95)&(df_merge['phastCons100way_mean']>0.2)
    mask4 = df_merge['GERP_max']>2
    mask5 = df_merge['AnnotSV_ACMG_class'].astype(float)>=4
    mask6 = df_merge['CADD-SV_PHRED-score']>10 
    mask7 = df_merge['AnnotSV_SV_type']=='TRA'
    mask8 = df_merge['Prob_random_segregation']<=(1/16)

    df_merge["prioritized_variant"]=0
    df_merge.loc[mask1|mask2|mask3|mask4|mask5|mask6|mask7|mask8,'prioritized_variant']=1
    df_merge.sort_values('prioritized_variant',ascending=False).to_csv(outcsv ,index=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--filt2_csv')
    parser.add_argument('--af')
    parser.add_argument('--vep')
    parser.add_argument('--cadd')
    parser.add_argument('--annotSV')
    parser.add_argument('--family_segregation')
    parser.add_argument('--master')
    parser.add_argument('--outcsv')
    args = parser.parse_args()
    main(args.filt2_csv, args.af, args.vep, args.cadd, args.annotSV, args.family_segregation, args.master, args.outcsv)

