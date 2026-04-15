snakemake all \
  --cores 8 \
  --dry-run \
  --resources mem_mb=15000 copyq=8 \
  --configfile config.yaml \
  --keep-going \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --use-conda \
  --printshellcmds 
