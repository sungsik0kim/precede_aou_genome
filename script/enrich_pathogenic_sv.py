import os
import pandas as pd
import gzip
from tqdm import tqdm
import numpy as np
import argparse


def main(csv, target, case, ctrl, outcsv):
    case_id = pd.read_csv(case)['analysis_sample_id'].astype(str).tolist()
    ctrl_id = pd.read_csv(ctrl)['analysis_sample_id'].astype(str).tolist()
    
    ###### Read Data #######
    target = pd.read_csv(target,header=None)[0].unique().tolist()
    df_var = pd.read_csv(csv)
    df_var.columns = [x.split('.sniffles2')[0] for x in df_var.columns.tolist()]
    id_ctrl= pd.read_csv(ctrl)["analysis_sample_id"].astype(str).tolist()
    id_case= pd.read_csv(case)["analysis_sample_id"].astype(str).tolist()
    
    ###### Parse Data #######
    df_var["END"] = df_var['INFO'].str.split(';').str[3].str.split('=').str[1]
    df_var.loc[df_var['INFO'].str.split(';').str[3].str.split('=').str[0]!='END','END']=''

    l_record = []
    for _, row in df_var.iterrows():
        tot_ctrl = np.isfinite(row[id_ctrl].astype(float).values).sum()
        N_ctrl = (row[id_ctrl].astype(float).values>=1).sum()
        ctrl_af = N_ctrl/tot_ctrl

        tot_case = np.isfinite(row[id_case].astype(float).values).sum()
        N_case = (row[id_case].astype(float).values>=1).sum()
        case_af = N_case/tot_case

        OR = ((N_case+0.01) * (tot_ctrl - N_ctrl+0.01)) / \
             ((N_ctrl+0.01) * (tot_case - N_case+0.01))


        l_record.append({
            'SYMBOL': row["SYMBOL"],
            'Conseqence':row['impact'],
            'SVTYPE':row['SVTYPE'],
            'SVLEN':row['SVLEN'],
            'CHROM':row['CHROM'],
            'POS':row['POS'],
            'END':row['END'],
            "FPC_sample_id": ', '.join(row[id_case].index[row[id_case]>=1]),
            "OR": OR,
            "FPC_AF": case_af,
            "AoU_AF": ctrl_af,
            "CoLoRdb_AF":row['CoLoRdb_AF'],
            "HPRC_AF":row['HPRC_AF'],
            "gnomAD_Max_PopMax_AF": row['gnomAD_Max_PopMax_AF'],
            "CADD-SV_PHRED": row['CADD-SV_PHRED-score'],
            "variant ID": row["ID"],
            "cluster_id": row['cluster_id']
        })
    df_result = pd.DataFrame(l_record)
    df_result["target_gene"]=df_result["SYMBOL"].apply(lambda x: x in target)
    
    ###### Filtering #######
    df_result = df_result.loc[(df_result['FPC_AF']>0)&(df_result['FPC_AF']>df_result['AoU_AF'])]
    print('Variants in FPC cohort:')
    print(df_result.shape)

    mask = ((df_result['gnomAD_Max_PopMax_AF']<0.01) & (df_result['CoLoRdb_AF']<0.01)& (df_result['HPRC_AF']<0.01))
    df_result = df_result.loc[mask]
    print('Filtering common variants in gnomAD, CoLoRdb, HPRC:')
    print(df_result.shape)

    df_result = df_result.loc[(df_result['AoU_AF']<0.01)]
    print('Filtering common variants in AoU:')
    print(df_result.shape)

    df_result = df_result.sort_values('OR',ascending=False)
    df_result.to_csv(outcsv, index=False)
  

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="enrich pathogenic variants")
    parser.add_argument('--csv')
    parser.add_argument('--target')
    parser.add_argument('--case')
    parser.add_argument('--ctrl')
    parser.add_argument('--outcsv')
    args = parser.parse_args()
    
    main(args.csv, args.target, args.case, args.ctrl, args.outcsv)
