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

def main(filt2_csv, dbnsfp, vep, promoterai, family_segregation, master, outcsv):
    ### Read Data ###
    df_filt2= pd.read_csv(filt2_csv)
    df_filt2_slim = df_filt2[[
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'Consequence',
        'IMPACT',
        'Gene',
        'Existing_variation',
        'gnomADg_Max_AF',
        'HPRC_AF',
        'CoLoRSdb_AF',
        'case_EAF',
        'ctrl_AF',
        'AoU_Max_AF',
        'internal_ctrl_AF',
        'Unaffected_Rel_EAF',
        'OR',
        'FPC_sample_id'    
    ]]
    df_filt2_slim = df_filt2_slim[~df_filt2_slim['FPC_sample_id'].isna()].reset_index(drop=True)
    
    df_dbnsfp = read_vcf(dbnsfp)
    df_vep = read_vcf(vep)
    df_promoterai = pd.read_csv(promoterai, sep="\t")
    df_segregation = pd.read_csv(family_segregation)
    df_master = pd.read_csv(master)

    #### Parse VEP ###
    df_vep["_vep"] = expand_col(df_vep['INFO'], 'CSQ')
    vep_cols='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|GIVEN_REF|USED_REF|BAM_EDIT|SOURCE|HGVS_OFFSET|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|SpliceAI_cutoff|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|CADD_PHRED|CADD_RAW|ClinVar|ClinVar_CLNSIG|ClinVar_CLNDN|phyloP100way|phyloP100way_min|phyloP100way_mean|phyloP100way_max|phastCons100way|phastCons100way_min|phastCons100way_mean|phastCons100way_max|GERP|GERP_min|GERP_mean|GERP_max'.split('|')
    for enum, col in enumerate(vep_cols):
        df_vep[col]= df_vep["_vep"].str.split('|').str[enum]
    vep_slim_cols = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'SYMBOL',
        'HGVSc',
        'HGVSp',
        'SpliceAI_cutoff',
        'SpliceAI_pred_DP_AG',
        'SpliceAI_pred_DP_AL',
        'SpliceAI_pred_DP_DG',
        'SpliceAI_pred_DP_DL',
        'SpliceAI_pred_DS_AG',
        'SpliceAI_pred_DS_AL',
        'SpliceAI_pred_DS_DG',
        'SpliceAI_pred_DS_DL',
        'SpliceAI_pred_SYMBOL',
        'ClinVar',
        'ClinVar_CLNSIG',
        'ClinVar_CLNDN',
        'phyloP100way_min',
        'phyloP100way_mean',
        'phyloP100way_max',
        'phastCons100way_min',
        'phastCons100way_mean',
        'phastCons100way_max',
        'GERP_min',
        'GERP_mean',
        'GERP_max',
        'CADD_PHRED',
        'CADD_RAW',
    ]
    df_vep_slim = df_vep[vep_slim_cols]
    cols = [
            'phyloP100way_min', 'phyloP100way_mean', 'phyloP100way_max',
            'phastCons100way_min','phastCons100way_mean', 'phastCons100way_max',
            'GERP_min', 'GERP_mean', 'GERP_max'                 
           ]
    for col in cols:
        df_vep_slim.loc[df_vep_slim[col]=='',col] = np.nan
        df_vep_slim[col] = df_vep_slim[col].astype(float)
    
    #### dbNSFP
    dbnsfp_tags='Interpro_domain,REVEL_score,REVEL_rankscore,VEST4_score,VEST4_rankscore,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,ESM1b_score,ESM1b_converted_rankscore,ESM1b_pred,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,VEP_canonical,MANE'.split(',')
    for tag in dbnsfp_tags:
        df_dbnsfp[tag] = expand_col(df_dbnsfp['INFO'], tag=tag)

    mask = ~df_dbnsfp['Ensembl_geneid'].isna()
    df_dbnsfp = df_dbnsfp.loc[mask]

    mask = mask & (df_dbnsfp['MANE'].apply(lambda x: "Select" in x) & (df_dbnsfp['VEP_canonical'].apply(lambda x: "YES" in x) ))
    df_dbnsfp = df_dbnsfp.loc[mask].copy().reset_index(drop=True)

    dbnsfp_tags='Interpro_domain,REVEL_score,VEST4_score,AlphaMissense_score,AlphaMissense_pred,ESM1b_score,ESM1b_pred,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,VEP_canonical,MANE'.split(',')
    l_record = []
    for _,row in df_dbnsfp.iterrows():
        mane = row['MANE'].split(',').index('Select')
        l_record.append({tag:row[tag].split(',')[mane] for tag in dbnsfp_tags})
    _df_dbnsfp = pd.DataFrame(l_record)
    for col in dbnsfp_tags:
        df_dbnsfp[col] = _df_dbnsfp[col]

    dbnsfp_slim_cols = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'Interpro_domain', 'REVEL_score', 'REVEL_rankscore', 'VEST4_score',
        'VEST4_rankscore', 'AlphaMissense_score', 'AlphaMissense_rankscore',
        'AlphaMissense_pred', 'ESM1b_score', 'ESM1b_converted_rankscore',
        'ESM1b_pred', 'Ensembl_geneid', 'Ensembl_transcriptid',
        'Ensembl_proteinid'
    ]
    df_dbnsfp_slim = df_dbnsfp[dbnsfp_slim_cols]
    
    ### PromoterAI ###
    df_promoterai = df_promoterai.rename(columns={'chrom':"CHROM"
                                              ,'pos':"POS"
                                              ,'ref':"REF"
                                              ,'alt':"ALT"
                                              ,'strand':"promoterAI_strand"
                                              ,'gene':"promoterAI_gene"
                                              ,'score':"PromoterAI_score"})
    
    ### Family Segregation ###
    df_segregation["CHROM"]=df_segregation["ID"].str.split('_').str[0]
    df_segregation["POS"]=df_segregation["ID"].str.split('_').str[1].astype(int)
    df_segregation["REF"]=df_segregation["ID"].str.split('_').str[2]
    df_segregation["ALT"]=df_segregation["ID"].str.split('_').str[3]
    df_segregation["Prob_random_segregation"]=df_segregation["N"]
    
    
    ### Merge ###
    df_merge = pd.merge(df_filt2_slim, df_vep_slim, on=['CHROM', 'POS', 'ID', 'REF', 'ALT'], how="inner")
    df_merge = pd.merge(df_merge, df_dbnsfp_slim, on=['CHROM', 'POS', 'ID', 'REF', 'ALT'], how="left")
    df_merge = pd.merge(df_merge, df_promoterai, on=['CHROM', 'POS', 'REF', 'ALT'], how="left")
    df_merge = pd.merge(df_merge, df_segregation.drop(columns=["ID","N"]), on=['CHROM', 'POS', 'REF', 'ALT'], how="left")
    
    df_merge.loc[df_merge["CADD_PHRED"]=='','CADD_PHRED']=np.nan
    df_merge["CADD_PHRED"] = df_merge["CADD_PHRED"].astype(float)
    df_merge.loc[df_merge["REVEL_score"]=='.','REVEL_score']=np.nan
    df_merge["REVEL_score"] = df_merge["REVEL_score"].astype(float)
    df_merge.loc[df_merge["VEST4_score"]=='.','VEST4_score']=np.nan
    df_merge["VEST4_score"] = df_merge["VEST4_score"].astype(float)
    df_merge.loc[df_merge["ESM1b_score"]=='.','ESM1b_score']=np.nan
    df_merge["ESM1b_score"] = df_merge["ESM1b_score"].astype(float)
    
    l_ambry_accession=[]
    for sample_ids in df_merge["FPC_sample_id"].astype(str):
        sample_ids = sample_ids.split(', ')
        _l_ambry_accession=[]
        for sample_id in sample_ids:
            print(sample_id)
            _l_ambry_accession.append(str(df_master.loc[df_master['Analysis_sample_id']==sample_id,'Ambry accession #'].iloc[0]))
        l_ambry_accession.append(', '.join(_l_ambry_accession))
    df_merge['Ambry_accession'] = l_ambry_accession
    df_merge['Ambry_accession'] = df_merge['Ambry_accession'].str.replace('nan','N/A')

    
    ### Prioritize & Print ###
    mask1 = df_merge['ClinVar_CLNSIG'].astype(str).apply(lambda x: ("pathogenic" in x) or ("Pathogenic" in x) )
    mask2 = df_merge["SpliceAI_cutoff"]=="PASS"
    mask3 = df_merge["Prob_random_segregation"]<=(1/16)
    mask4 = df_merge["REVEL_score"]>0.8
    mask5 = df_merge["VEST4_score"]>0.8
    mask6 = df_merge['AlphaMissense_pred']=="P"
    mask7 = df_merge['ESM1b_score']< -10
    mask8 = df_merge['PromoterAI_score'].abs()>0.5
    mask9 = df_merge['phyloP100way_max']>7
    mask10 = df_merge['GERP_max']>3


    df_merge["prioritized_variant"]=0
    df_merge.loc[mask1|mask2|mask3|mask4|mask5|mask6|mask7|mask8|mask9|mask10,'prioritized_variant']=1
    df_merge.sort_values('prioritized_variant',ascending=False).to_csv(outcsv,index=False)

    
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--filt2_csv')
    parser.add_argument('--dbnsfp')
    parser.add_argument('--vep')
    parser.add_argument('--promoterai')
    parser.add_argument('--family_segregation')
    parser.add_argument('--master')
    parser.add_argument('--outcsv')
    args = parser.parse_args()
    main(args.filt2_csv, args.dbnsfp, args.vep, args.promoterai, args.family_segregation, args.master, args.outcsv)

