snakemake all \
  --cores 16 \
  --resources mem_mb=52000 copyq=8 \
  --configfile config.yaml \
  --keep-going \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --use-conda \
  --printshellcmds 
