#!/bin/bash
# Usage: run_vep.sh input_vcf output_dir output_prefix

VCF_FILE="$1" # should be sorted and indexed (.vcf.gz)
OUTDIR="$2"
PREFIX="$3"
PLUGIN="/home/jupyter/.vep/Plugins"
CLINVAR="/home/jupyter/workspaces/longreadseqcontrolset/reference/ClinVar"
PHYLOP="/home/jupyter/workspaces/longreadseqcontrolset/reference/phyloP"
PHASTCONS="/home/jupyter/workspaces/longreadseqcontrolset/reference/phastCons"
GERP="/home/jupyter/workspaces/longreadseqcontrolset/reference/GERP"


# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

#### VEP ####
/home/jupyter/workspaces/longreadseqcontrolset/bin/ensembl-vep/vep --dir /home/jupyter/.vep \
  --fasta /home/jupyter/workspaces/longreadseqcontrolset/reference/human_GRCh38_no_alt_analysis_set.fasta \
  --assembly GRCh38 --cache --offline --no_stats --format vcf --vcf --force_overwrite --max_sv_size 10000000 --distance 10000 \
  --fork 16 --buffer_size 200000 \
  --regulatory --check_existing --pubmed --hgvs \
  --symbol --canonical --pick \
  --custom file="$CLINVAR/clinvar.vcf.gz",short_name=ClinVar,format=vcf,type=overlap,coords=0,fields=CLNSIG%CLNDN \
  --custom file="$PHYLOP/hg38.phyloP100way.bw",short_name=phyloP100way,format=bigwig,type=overlap,coords=0,summary_stats=min%mean%max \
  --custom file="$PHASTCONS/hg38.phastCons100way.bw",short_name=phastCons100way,format=bigwig,type=overlap,coords=0,summary_stats=min%mean%max \
  --custom file="$GERP/gerp_conservation_scores.homo_sapiens.GRCh38.bw",short_name=GERP,format=bigwig,type=overlap,coords=0,summary_stats=min%mean%max \
  --input_file "$VCF_FILE" \
  --output_file "$OUTDIR/$PREFIX.vep.vcf"

  

