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

def calculate_allele_counts(df, cols):
    """
    Logically robust allele counting.
    Counts '1's specifically, handling phased (|) and unphased (/) 
    and multi-allelics (1/2) correctly for Allele 1.
    """
    # 1. Get the relevant dataframe subset
    dsub = df[cols]
    
    # 2. Calculate Denominator (Total Alleles)
    # Count samples that are NOT missing (./. or .)
    # Note: simple inequality checks can be risky with complex VCF strings, 
    # but strictly speaking, ./. is the standard missing GT.
    called_samples = (dsub != './.') & (dsub != '.')
    # Total alleles = count of called samples * 2 (assuming diploid)
    n_alleles = called_samples.sum(axis=1) * 2
    
    # 3. Calculate Numerator (Alt Alleles)
    # We count how many times '1' appears in the genotype string.
    # '0/1' -> 1
    # '1/1' -> 2
    # '1/2' -> 1 (Correctly counts Allele 1, ignores Allele 2)
    # '0/2' -> 0 (Correctly ignores Allele 2)
    # ./1   -> 1 (Handles half-calls if they exist)
    n_alt = dsub.apply(lambda x: x.str.count('1')).sum(axis=1)
    return n_alleles, n_alt


def expand_col(series, tag, dtype=None):
    if dtype:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0].astype(dtype)
    else:
        return series.str.split(f'{tag}=').str[1].str.split(';').str[0]

def get_EAF (df_genotype, sample_id, pedigree_id):
    df_info = pd.DataFrame({
        'sample_id':sample_id
        ,'pedigree_id':pedigree_id
    })
    df_info_singleton = df_info.loc[df_info['pedigree_id'].isna()]
    df_info_non_single = df_info.loc[~df_info['pedigree_id'].isna()]

    ## For Singletons
    cols = [col for col in df_info_singleton['sample_id'].tolist() if col in df_genotype.columns]
    base, alt = calculate_allele_counts(df_genotype, cols)

    for _, df_gb in df_info_non_single.groupby('pedigree_id'):
        cols = [col for col in df_gb['sample_id'].tolist() if col in df_genotype.columns]
        if cols:
            _base, _alt = calculate_allele_counts(df_genotype, cols)
            _base = _base/len(cols)
            _alt = _alt/len(cols)

            base+=_base
            alt+=_alt

    af = alt/base
    return af

def main(vcf, af, case, ctrl, ctrl2, ctrl3, master, aou_covar, outdir):
    print('###### Read and Parse Data #######')
    case_id = pd.read_csv(case)['analysis_sample_id'].astype(str).tolist()
    ctrl_id = pd.read_csv(ctrl)['analysis_sample_id'].astype(str).tolist()   

    df_af = pd.read_csv(af)
    
    ## VCF for Precede FPC and AoU
    df_vcf = read_vcf_gz(vcf)
    df_int_ctrl = read_vcf_gz(ctrl2)
    df_unaffected = read_vcf_gz(ctrl3)
    df_master = pd.read_csv(master)
    df_aou_covar = pd.read_csv(aou_covar,sep="\t")

    print('## Step 1: Filtering by AF')
    ## Filtering by Databases ##
    mask_gnomad = (df_af['gnomAD_Max_PopMax_AF']<0.01)
    mask_colors = (df_af['CoLoRdb_AF']<0.01)
    mask_hprc = (df_af['HPRC_AF']<0.01)

    ## Filtering by AoU + Precede ##
    case_id2 = [col for col in case_id if col in df_vcf.columns.tolist()]
    ctrl_id2 = [f'{col}.sniffles2.hg38' for col in ctrl_id if f'{col}.sniffles2.hg38' in df_vcf.columns.tolist()]
    for col in tqdm(case_id2+ctrl_id2):
        df_vcf[col] = df_vcf[col].str.split(':').str[0]
    
    ## -- EAF for case
    df_vcf['case_EAF'] = get_EAF(df_vcf[case_id2]
                                 ,df_master['Analysis_sample_id'].tolist()
                                 ,df_master['Analysis_pedigree_id'].tolist()
                                )
    ## -- AF for AoU Ctrl
    base, alt = calculate_allele_counts(df_vcf, ctrl_id2)
    af = alt/base
    df_vcf['ctrl_AF'] = af
    
    ## Filtering by Precede Internal Ctrl ##
    sample_cols = df_int_ctrl.columns[9:]
    for col in sample_cols:
        df_int_ctrl[col] = df_int_ctrl[col].str.split(':').str[0]
    base, alt = calculate_allele_counts(df_int_ctrl, sample_cols)
    af = alt/base
    df_int_ctrl['internal_ctrl_AF'] = af
    df_vcf = pd.merge(df_vcf, df_int_ctrl[['ID','internal_ctrl_AF']], on = 'ID', how='left')
    
    ## Filtering by Precede Unaffected Rel (only Case cohort independent samples used) ##
    case_pedigree = df_master.loc[df_master['Analysis_sample_id'].apply(lambda x: x in case_id),'Analysis_pedigree_id'].unique()
    case_unrelated_ids = df_master.loc[df_master['Analysis_pedigree_id'].apply(lambda x: x not in case_pedigree), 'Analysis_sample_id'].tolist()
    case_unrelated_ids = [id for id in case_unrelated_ids if id not in case_id] 
    sample_cols = df_unaffected.columns[9:]
    sample_cols = [col for col in sample_cols if col in case_unrelated_ids]
    for col in [col for col in sample_cols]:
        df_unaffected[col] = df_unaffected[col].str.split(':').str[0]
    df_unaffected['Unaffected_Rel_EAF'] = get_EAF(df_genotype = df_unaffected[sample_cols]
                                            , sample_id = df_master['Analysis_sample_id'].tolist()
                                            , pedigree_id =df_master['Analysis_pedigree_id'].tolist()
                                           )
    df_vcf = pd.merge(df_vcf, df_unaffected[['ID','Unaffected_Rel_EAF']], on = 'ID', how='left')

    ## Filtering 
    ## ::Rare in AoU?
    l_eur = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='eur', 'person_id'].tolist()
    l_afr = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='afr', 'person_id'].tolist()
    l_amr = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='amr', 'person_id'].tolist()
    l_eas = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='eas', 'person_id'].tolist()
    l_sas = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='sas', 'person_id'].tolist()
    l_mid = df_aou_covar.loc[df_aou_covar['genetic_ancestry']=='mid', 'person_id'].tolist()

    df_vcf_aou = df_vcf[ctrl_id2].copy()
    df_vcf_aou.columns = [int(col.split('.sniffles2')[0]) for col in df_vcf[ctrl_id2].columns]
    df_vcf_aou_eur = df_vcf_aou[[s for s in l_eur if s in df_vcf_aou.columns]]
    df_vcf_aou_afr = df_vcf_aou[[s for s in l_afr if s in df_vcf_aou.columns]]
    df_vcf_aou_amr = df_vcf_aou[[s for s in l_amr if s in df_vcf_aou.columns]]
    df_vcf_aou_eas = df_vcf_aou[[s for s in l_eas if s in df_vcf_aou.columns]]
    df_vcf_aou_sas = df_vcf_aou[[s for s in l_sas if s in df_vcf_aou.columns]]
    df_vcf_aou_mid = df_vcf_aou[[s for s in l_mid if s in df_vcf_aou.columns]]

    df_vcf.loc[:,'AoU_EUR_AF'] = get_EAF(df_vcf_aou_eur, df_vcf_aou_eur.columns, df_vcf_aou_eur.columns)
    df_vcf.loc[:,'AoU_AFR_AF'] = get_EAF(df_vcf_aou_afr, df_vcf_aou_afr.columns, df_vcf_aou_afr.columns)
    df_vcf.loc[:,'AoU_AMR_AF'] = get_EAF(df_vcf_aou_amr, df_vcf_aou_amr.columns, df_vcf_aou_amr.columns)
    df_vcf.loc[:,'AoU_EAS_AF'] = get_EAF(df_vcf_aou_eas, df_vcf_aou_eas.columns, df_vcf_aou_eas.columns)
    df_vcf.loc[:,'AoU_SAS_AF'] = get_EAF(df_vcf_aou_sas, df_vcf_aou_sas.columns, df_vcf_aou_sas.columns)
    df_vcf.loc[:,'AoU_MID_AF'] = get_EAF(df_vcf_aou_mid, df_vcf_aou_mid.columns, df_vcf_aou_mid.columns)
    df_vcf.loc[:,'AoU_Max_AF'] = df_vcf[['AoU_EUR_AF','AoU_AFR_AF','AoU_AMR_AF','AoU_EAS_AF','AoU_SAS_AF','AoU_MID_AF']].max(axis=1)

    mask_aou_1 = (df_vcf['AoU_Max_AF']<0.01)
#     mask_aou_2 = (df_vcf['case_EAF']>(df_vcf['AoU_Max_AF'])) # AF@Case > AF@AoU

    ## ::Higher AF than Non-FPC samples? :: Based on Effective Allele Frequenct (EAF)
    mask_precede_1 = (df_vcf['case_EAF']>(df_vcf['internal_ctrl_AF'])) # AF@Case > AF@Internal Ctrl
#     mask_precede_2 = (df_vcf['case_EAF']>(df_vcf['Unaffected_Rel_EAF'])) # AF@Case > AF@Unaffected Rel
    

    print('## Logging')
    print("Original Num of Variants : {}".format(df_vcf.shape[0]))
    print("gnomAD < 1% : {}".format(df_vcf.loc[mask_gnomad].shape[0]))
    print("HPRC<1% & CoLoRS<1% : {}".format(df_vcf.loc[mask_gnomad&mask_colors&mask_hprc].shape[0]))
    print("AoU<1% : {}".format(df_vcf.loc[mask_gnomad&mask_colors&mask_hprc&mask_aou_1].shape[0]))
    mask_db = mask_gnomad&mask_colors&mask_hprc&mask_aou_1
    print("FPC > internal_ctrl_AF : {}".format(df_vcf.loc[mask_db&mask_precede_1].shape[0]))
    
#     print("FPC > AoU : {}".format(df_vcf.loc[mask_db&mask_aou_1&mask_aou_2].shape[0]))
#     print("FPC > Unaffected_Rel_EAF : {}".format(df_vcf.loc[mask_db&mask_aou&mask_precede].shape[0]))
    
    
    print('###### Add columns: OR, FPC_sample_id #######')
    df_vcf['OR'] = (1/(1/df_vcf['case_EAF']-1))/(1/(1/df_vcf['ctrl_AF']-1))
    l_case = []
    for _, row in tqdm(df_vcf.iterrows()):
        gt=row[case_id]
        mask = ((gt!='0/0') & (gt!='0|0') & (gt!='./.'))
        l_case.append(', '.join(gt[mask].index.tolist()))
    df_vcf['FPC_sample_id'] = l_case
    df_filt = df_vcf.loc[mask_db&mask_precede_1].reset_index(drop=True) ## with filter
#     df_filt = df_vcf.reset_index(drop=True) ## no filter
    
    print('###### Save data ######')
    df_vcf.to_csv(os.path.join(outdir, "sv_filter2_before_filtering.csv"), index=False)
    df_filt.to_csv(os.path.join(outdir, "sv_filter2.csv"), index=False)
    df_filt[['ID']].to_csv(os.path.join(outdir, "sv_filter2.id"), index=False)
    df_filt_bed=df_filt[["CHROM","POS"]].copy()
    df_filt_bed["START"] = df_filt_bed["POS"]-1
    df_filt_bed["END"] = df_filt_bed["POS"]
    df_filt_bed = df_filt_bed[["CHROM","START","END"]]
    df_filt_bed.to_csv(os.path.join(outdir, "sv_filter2.bed"), index=False, sep="\t",header=None)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="enrich pathogenic variants")
    parser.add_argument('--vcf')
    parser.add_argument('--af')
    parser.add_argument('--case')
    parser.add_argument('--ctrl')
    parser.add_argument('--ctrl2')
    parser.add_argument('--ctrl3')
    parser.add_argument('--master')
    parser.add_argument('--aou_covar')
    parser.add_argument('--outdir')
    args = parser.parse_args()
    main(args.vcf, args.af, args.case, args.ctrl, args.ctrl2, args.ctrl3, args.master, args.aou_covar, args.outdir)
