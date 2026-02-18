import argparse
import pandas as pd
from pathlib import Path
import os

def main(aou_snf_dir, precede_snf_dir, precede_mastersheet, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    ctrl_sample_out = os.path.join(outdir,'ctrl_sample.csv')
    case_sample_out = os.path.join(outdir,'case_sample.csv')
    internal_ctrl_sample_out= os.path.join(outdir,'internal_ctrl_sample.csv')
    unaffected_rel_sample_out= os.path.join(outdir,'unaffected_rel_sample.csv')
    
    # FPC cases
    df_master = pd.read_csv(precede_mastersheet)
    mask1 = df_master['dataset']!='normal ctrl'
    mask2 = df_master['Exclusion'].isna()
    mask3 = df_master['Affected Status']=="Affected"
    fpc = df_master.loc[mask1&mask2&mask3,'Analysis_sample_id'].astype(str).tolist()
    
    # internal controls
    mask1 = df_master['dataset']=='normal ctrl'
    mask2 = df_master['Exclusion'].isna()
    i_ctrl = df_master.loc[mask1&mask2,'Analysis_sample_id'].astype(str).tolist()
    
    # Unaffected family members
    mask1 = df_master['dataset']!='normal ctrl'
    mask2 = df_master['Exclusion'].isna()
    mask3 = df_master['Affected Status']=="Unaffected"
    unaffected = df_master.loc[mask1&mask2&mask3,'Analysis_sample_id'].astype(str).tolist()
    

    # For AoU
    l_ctrl_snf = Path(aou_snf_dir).glob('*.snf')
    l_sample = []
    l_path = []
    for path in l_ctrl_snf:
        l_path.append(path)
        l_sample.append(path.stem.split('.')[0])
    pd.DataFrame({
        'analysis_sample_id':l_sample,
        'path':l_path
    }).to_csv(ctrl_sample_out, index=False)    
    
    # For Precede
    l_precede_snf = list(Path(precede_snf_dir).glob("*.snf"))
    
    ## -- FPC
    l_sample = []
    l_path = []
    for path in l_precede_snf:
        analysis_sample_id = path.stem
        if analysis_sample_id in fpc:
            l_path.append(path)
            l_sample.append(analysis_sample_id)
    pd.DataFrame({
        'analysis_sample_id':l_sample,
        'path':l_path
    }).to_csv(case_sample_out, index=False)
    
    ## -- internal control
    l_sample = []
    l_path = []
    for path in l_precede_snf:
        analysis_sample_id = path.stem
        if analysis_sample_id in i_ctrl:
            l_path.append(path)
            l_sample.append(analysis_sample_id)
    pd.DataFrame({
        'analysis_sample_id':l_sample,
        'path':l_path
    }).to_csv(internal_ctrl_sample_out, index=False)
    
    ## -- Unaffected member
    l_sample = []
    l_path = []
    for path in l_precede_snf:
        analysis_sample_id = path.stem
        if analysis_sample_id in unaffected:
            l_path.append(path)
            l_sample.append(analysis_sample_id)
    pd.DataFrame({
        'analysis_sample_id':l_sample,
        'path':l_path
    }).to_csv(unaffected_rel_sample_out, index=False)
    
#         python {params.script} --aou_snf {input.aou_snf_dir} --precede_snf {input.precede_snf_dir} --precede_master {input.precede_mastersheet} --outdir {params.outdir} > {log} 2>&1


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Merge and Filter Variants")
    parser.add_argument('--aou_snf')
    parser.add_argument('--precede_snf')
    parser.add_argument('--precede_master')
    parser.add_argument('--outdir')
    args = parser.parse_args()
    
    main(args.aou_snf, args.precede_snf, args.precede_master, args.outdir)
    
