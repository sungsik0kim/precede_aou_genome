import argparse
import pandas as pd
from pathlib import Path
import os

def main(target_ld, target_gene, transcript, output):
    df_ld = pd.read_csv(target_ld, sep="\t")
    df_gene = pd.read_csv(target_gene, header=None)
    df_transcript =  pd.read_csv(transcript, sep="\t", header=None)
    
    s = df_transcript[8].str.split(';')
    df_transcript["gene_id"] = s.str[0].str.strip().str.split().str[1].str[1:-1]
    df_transcript["gene_name"] = s.str[3].str.strip().str.split().str[1].str[1:-1]
    df_transcript = df_transcript[["gene_name","gene_id",0,3,4,6]]
    df_transcript.columns = ['gene_name','gene_id','chr','start','end','strand']
    df_transcript["chrom"] = df_transcript["chr"].str[3:]
    
    ## 200kb covers common short-range enhancer-promoter interactions
    df_temp = pd.merge(df_transcript, df_ld, on='chrom')
    mask1 = df_temp["start"] < df_temp["block_end"]+200000
    mask2 = df_temp["end"] > df_temp["block_start"]-200000
    l_genes = df_temp[mask1&mask2]["gene_name"].unique().tolist()
    l_genes.extend(df_gene[0].tolist())
    l_genes = list(set(l_genes))
    
    with open(output, 'w') as f:
        for g in l_genes: f.write(g+'\n')
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="make gene list combining LD_blocks")
    parser.add_argument('--target_ld')
    parser.add_argument('--target_gene')
    parser.add_argument('--transcript')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(args.target_ld, args.target_gene, args.transcript , args.output)
    
