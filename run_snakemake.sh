snakemake all \
  --dry-run \
  --cores 32 \
  --resources mem_mb=104000 copyq=8 \
  --configfile config.yaml \
  --keep-going \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --use-conda \
  --printshellcmds 
