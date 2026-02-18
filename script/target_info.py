import argparse
import pandas as pd
from pathlib import Path
import os

def main(target, transcript, outdir, padding=10000):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    t = pd.read_csv(target,header=None)[0].tolist()
    
    ## Transcript
    df = pd.read_csv(transcript, sep="\t", header=None)
    s = df[8].str.split(';')
    df["gene_id"] = s.str[0].str.strip().str.split().str[1].str[1:-1]
    df["gene_name"] = s.str[3].str.strip().str.split().str[1].str[1:-1]
    df = df[["gene_name","gene_id",0,3,4,6]]
    df.columns = ['gene_name','gene_id','chr','start','end','strand']
    df= df.loc[df['gene_name'].apply(lambda x: x in t)]
    missing = set(t) - set(df["gene_name"])
    if missing: raise ValueError(f"Gene info failed to search: {sorted(missing)}")
    df.to_csv(os.path.join(outdir,'transcript.csv'),index=False)
    
    df_target = df[["chr","start","end","gene_name"]].copy()
    df_target["start"]=df_target["start"]-padding
    df_target.loc[df_target["start"]<0 ,"start"]=0
    df_target["end"]=df_target["end"]+padding
    df_target.to_csv(os.path.join(outdir,'target_gene.bed'), index=False, header=False, sep="\t")
    
    # Transcript to Promoter
    padding_promoter = 10000
    l_record = []
    for _, row in df.iterrows():
        chr= row['chr']
        gene_name = row['gene_name']
        _start = row['start']
        _end = row['end']
        strand = row['strand']
        if strand =='+':
            start = max(0,_start-padding_promoter)
            end = _start+2000
            l_record.append({
                'chr':chr
                ,'start':start
                ,'end':end
                ,'gene_name':gene_name
                ,'strand':strand
            })
        if strand =='-':
            start = max(0,_end-2000)
            end = _end+padding_promoter
            l_record.append({
                'chr':chr
                ,'start':start
                ,'end':end
                ,'gene_name':gene_name
                ,'strand':strand
            })

    df_promoter = pd.DataFrame(l_record)
    df_promoter.to_csv(os.path.join(outdir,'target_gene_promoter.bed'), index=False, header=False, sep="\t")
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Merge and Filter Variants")
    parser.add_argument('--target')
    parser.add_argument('--transcript')
    parser.add_argument('--outdir')
    parser.add_argument('--padding', default=10000, type=int)
    args = parser.parse_args()
    main(args.target, args.transcript, args.outdir, args.padding)
    
