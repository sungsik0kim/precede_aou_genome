#!/bin/bash
# Usage: run_vep.sh input_vcf output_dir output_prefix

VCF_FILE="$1" # should be sorted and indexed (.vcf.gz)
OUTDIR="$2"
PREFIX="$3"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

#### VEP ####
/home/jupyter/workspaces/longreadseqcontrolset/bin/ensembl-vep/vep --dir /home/jupyter/.vep \
  --fasta /home/jupyter/workspaces/longreadseqcontrolset/reference/human_GRCh38_no_alt_analysis_set.fasta \
  --assembly GRCh38 --cache --offline --no_stats --format vcf --vcf --force_overwrite --pick --dist 100 \
  --fork 16 --buffer_size 200000 \
  --af_gnomadg \
  --input_file "$VCF_FILE" \
  --output_file "$OUTDIR/$PREFIX.gnomad.vcf"
