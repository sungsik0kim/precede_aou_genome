#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import pysam

REF_FA = '/home/jupyter/workspaces/longreadseqcontrolset/reference/human_GRCh38_no_alt_analysis_set.fasta'

def get_sequence_from_ref(chrom, pos_1based, indel_length):
    with pysam.FastaFile(REF_FA) as ref_genome:
        # Convert 1-based pos to 0-based for pysam
        start_0based = pos_1based - 1 
        
        # Fetch bases: 1 anchor base + indel_length
        ref_seq_full = ref_genome.fetch(chrom, start_0based, start_0based + indel_length + 1)
    return ref_seq_full

def process_sv(input_path, output_path):
    print(f"Loading input file: {input_path}")
    df_sv = pd.read_csv(input_path, low_memory=False)

    # Filter to limit length < 200k
    df_sv = df_sv[df_sv['AnnotSV_SV_length'].abs() < 200000] 
    
    cols = ["CHROM", "POS", "REF", "ALT", "AnnotSV_SV_length"]
    df_sv_ins = df_sv.loc[df_sv['AnnotSV_SV_type'] == "INS", cols]
    df_sv_del = df_sv.loc[df_sv['AnnotSV_SV_type'] == "DEL", cols]
    df_sv_inv = df_sv.loc[df_sv['AnnotSV_SV_type'] == "INV", cols]
    df_sv_dup = df_sv.loc[df_sv['AnnotSV_SV_type'] == "DUP", cols]

    print("Processing Insertions...")
    df_sv_ins["REF"] = [get_sequence_from_ref(chrom, pos, 0) for chrom, pos in zip(df_sv_ins["CHROM"].tolist(), df_sv_ins["POS"].astype(int).tolist())]
    df_sv_ins["ALT"] = df_sv_ins["REF"] + df_sv_ins["ALT"].astype(str)
    
    # Random sequence assignment to symbolic insertions
    idx_symbolic = df_sv_ins[df_sv_ins["ALT"].apply(lambda x: "<INS>" in str(x))].index
    for idx in idx_symbolic:
        n = df_sv_ins.loc[idx, "AnnotSV_SV_length"].astype(int)
        df_sv_ins.loc[idx, "ALT"] = ''.join(np.random.choice(["A", "G", "C", "T"], n))

    print("Processing Deletions...")
    sequences = [get_sequence_from_ref(chrom, pos, svlen) for chrom, pos, svlen \
                        in zip(df_sv_del["CHROM"].tolist()
                               , df_sv_del["POS"].astype(int).values - 1
                               , df_sv_del["AnnotSV_SV_length"].abs().astype(int).values
                              )
                ]
    df_sv_del["REF"] = [seq for seq in sequences]
    df_sv_del["ALT"] = [seq[0] for seq in sequences]

    print("Processing Inversions...")
    sequences = [get_sequence_from_ref(chrom, pos, svlen) for chrom, pos, svlen \
                        in zip(df_sv_inv["CHROM"].tolist()
                               , df_sv_inv["POS"].astype(int).values
                               , df_sv_inv["AnnotSV_SV_length"].abs().astype(int).values - 1
                              )
                ]
    trans_table = str.maketrans("ATGC", "TACG")
    df_sv_inv["REF"] = sequences
    df_sv_inv["ALT"] = [seq.translate(trans_table)[::-1] for seq in sequences]

    print("Processing Duplications...")
    sequences = [get_sequence_from_ref(chrom, pos, svlen) for chrom, pos, svlen \
                        in zip(df_sv_dup["CHROM"].tolist()
                               , df_sv_dup["POS"].astype(int).values
                               , df_sv_dup["AnnotSV_SV_length"].abs().astype(int).values - 1
                              )
                ]
    df_sv_dup["REF"] = sequences
    df_sv_dup["ALT"] = [seq + seq for seq in sequences]

    print("Concatenating into df_sv_slim...")
    df_sv_slim = pd.concat((df_sv_ins, df_sv_del, df_sv_inv, df_sv_dup))[["CHROM", "POS", "REF", "ALT"]].reset_index(drop=True)
    
    # Add variant_str feature (from cell 6 of your notebook)
    df_sv_slim['variant_str'] = df_sv_slim['CHROM'].astype(str) + ':' + df_sv_slim['POS'].astype(str) + ':' + df_sv_slim['REF'].astype(str) + ':' + df_sv_slim['ALT'].astype(str)
    df_sv_slim["dummy"]=""

    # Save to CSV
    df_sv_slim[["variant_str","dummy"]].to_csv(output_path, index=False, header=None)
    print(f"Success! Output saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform sv_prioritize.csv to a AG variant format.")
    parser.add_argument("-i", "--input", default="input/sv_prioritize.csv", help="Path to input sv_prioritize.csv")
    parser.add_argument("-o", "--output", default="output/df_sv_slim.csv", help="Path for the output CSV file")
    
    args = parser.parse_args()
    process_sv(args.input, args.output)