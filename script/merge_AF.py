import argparse
import pandas as pd

def main(gnomAD, CoLoRdb, hprc, output):
    
    ### Read AF annotated vcf files
    skip_lines = 0
    with open(gnomAD, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_gnomad = pd.read_csv(
        gnomAD, 
        sep='\t', 
        skiprows=skip_lines
    )
    df_gnomad = df_gnomad.rename(columns={'#CHROM': 'CHROM'})

    skip_lines = 0
    with open(CoLoRdb, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_colordb = pd.read_csv(
        CoLoRdb, 
        sep='\t', 
        skiprows=skip_lines
    )
    df_colordb = df_colordb.rename(columns={'#CHROM': 'CHROM'})

    skip_lines = 0
    with open(hprc, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                skip_lines += 1
            else:
                break
    df_hprc = pd.read_csv(
        hprc, 
        sep='\t', 
        skiprows=skip_lines
    )
    df_hprc = df_hprc.rename(columns={'#CHROM': 'CHROM'})
    
    ### Parsing AFs
    df_gnomad["gnomAD_Max_PopMax_AF"] = df_gnomad["INFO"].str.split('Max_PopMax_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_AFR_AF"] = df_gnomad["INFO"].str.split('Max_AFR_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_AMI_AF"] = df_gnomad["INFO"].str.split('Max_AMI_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_AMR_AF"] = df_gnomad["INFO"].str.split('Max_AMR_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_ASJ_AF"] = df_gnomad["INFO"].str.split('Max_ASJ_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_EAS_AF"] = df_gnomad["INFO"].str.split('Max_EAS_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_EUR_AF"] = df_gnomad["INFO"].str.split('Max_EUR_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_FIN_AF"] = df_gnomad["INFO"].str.split('Max_FIN_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_MID_AF"] = df_gnomad["INFO"].str.split('Max_MID_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_NFE_AF"] = df_gnomad["INFO"].str.split('Max_NFE_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_OTH_AF"] = df_gnomad["INFO"].str.split('Max_OTH_AF=').str[1].str.split(';').str[0].astype(float)
    df_gnomad["gnomAD_Max_SAS_AF"] = df_gnomad["INFO"].str.split('Max_SAS_AF=').str[1].str.split(';').str[0].astype(float)

    df_colordb['CoLoRdb_AF'] = df_colordb["INFO"].str.split('CoLoRdb_AF=').str[1].str.split(';').str[0].astype(float)
    
    df_hprc['HPRC_AF'] = df_hprc["INFO"].str.split('hprc_AF=').str[1].str.split(';').str[0].astype(float)
    
    ### Merging AFs
    gnomad_cols = ['CHROM', 'POS', 'ID'
                   ,'gnomAD_Max_PopMax_AF', 'gnomAD_Max_AFR_AF', 'gnomAD_Max_AMI_AF', 'gnomAD_Max_AMR_AF'
                   , 'gnomAD_Max_ASJ_AF', 'gnomAD_Max_EAS_AF', 'gnomAD_Max_EUR_AF', 'gnomAD_Max_FIN_AF'
                   , 'gnomAD_Max_MID_AF', 'gnomAD_Max_NFE_AF', 'gnomAD_Max_OTH_AF', 'gnomAD_Max_SAS_AF']
    colordb_cols = ['CHROM', 'POS', 'ID','CoLoRdb_AF']
    hprc_cols = ['CHROM', 'POS', 'ID','HPRC_AF']
    df_merge = pd.merge(df_gnomad[gnomad_cols], df_colordb[colordb_cols], on = ['CHROM', 'POS', 'ID'], how='inner')
    df_merge = pd.merge(df_merge, df_hprc[hprc_cols], on = ['CHROM', 'POS', 'ID'], how='inner')
    df_merge = df_merge.fillna(0)
    
    df_merge.to_csv(output, index=False)
    
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="python merge_AF.py --gnomAD --CoLoRdb --hprc --output")
    parser.add_argument('--gnomAD')
    parser.add_argument('--CoLoRdb')
    parser.add_argument('--hprc')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(args.gnomAD, args.CoLoRdb, args.hprc, args.output)


