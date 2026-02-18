import argparse
import pandas as pd
import gzip

def main(vcf, promoter_bed, output):
    skip_lines = 0
    with gzip.open(vcf, 'rt') as f:
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
    
    df_bed = pd.read_csv(
        promoter_bed, 
        sep='\t', 
        header=None
    )

    df_out = df_vcf[['#CHROM','POS','REF','ALT']].copy().rename(columns={
        '#CHROM':'chrom'
        ,'POS':'pos'
        ,'REF':'ref'
        ,'ALT':'alt'
    })
    
    l_record =[] 
    for _, row in df_out.iterrows():
        mask = (row['chrom']==df_bed[0])
        mask = mask & (row['pos'] >= df_bed[1]) & (row['pos'] <= df_bed[2])
        for idx in df_bed.loc[mask].index:
            row_promoter = df_bed.loc[idx]
            l_record.append({
                'chrom': row['chrom']
                ,'pos': row['pos']
                ,'ref': row['ref']
                ,'alt': row['alt']
                ,'strand': row_promoter[4]
                ,'gene': row_promoter[3]
            })
    
    pd.DataFrame(l_record).to_csv(output, sep="\t", index=False)
    
if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Prep promoterAI input tsv")
    parser.add_argument('--vcf')
    parser.add_argument('--promoter_bed')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(args.vcf, args.promoter_bed, args.output)

