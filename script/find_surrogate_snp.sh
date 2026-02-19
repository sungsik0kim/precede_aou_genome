#!/bin/bash
# Usage: find_surrogate_snp.sh {sv_dir} {snv_dir} {outdir}

SV_DIR=$1
SNV_DIR=$2
OUTDIR=$3
RENAME=$4

mkdir -p $OUTDIR/tmp
################################
echo "PREP FOR SNV/INDEL"
################################

echo "Make BED region"
Script='/home/jupyter/workspaces/longreadseqcontrolset/git/precede_aou_tier1/script/target_info.py'
python $Script --target /home/jupyter/workspaces/longreadseqcontrolset/git/precede_aou_tier1/resource/target_gene.txt \
               --transcript /home/jupyter/workspaces/longreadseqcontrolset/output/precede_aou_tier1/target_info/transcript.gtf \
               --outdir $OUTDIR \
               --padding 500000

sort -k1,1V -k2,2n $OUTDIR/target_gene.bed > $OUTDIR/target_gene.sorted.bed
bgzip -c $OUTDIR/target_gene.sorted.bed > $OUTDIR/target_gene.sorted.bed.gz
tabix -p bed $OUTDIR/target_gene.sorted.bed.gz

echo "Intersect the SNV/Indels using the BED"
for CHROM in {1..22} X Y; do
    bcftools view --threads 4 -R $OUTDIR/target_gene.sorted.bed.gz -Oz \
        -o "$OUTDIR/tmp/snv_intersect_chr${CHROM}.vcf.gz" \
        $SNV_DIR/AoU_Precede_chr"${CHROM}".vcf.gz &
done
wait

echo "Index Intersected VCFs"
for chr in {1..22} X Y; do tabix -p vcf "$OUTDIR/tmp/snv_intersect_chr${chr}.vcf.gz" & done; wait


echo "Combine Chromosomes"
bcftools concat $OUTDIR/tmp/snv_intersect_chr*.vcf.gz -Oz -o "$OUTDIR/tmp/snv_intersect_unsorted.vcf.gz" --thread 8
bcftools sort -Oz -o $OUTDIR/snv.vcf.gz --temp-dir "$OUTDIR/tmp/" "$OUTDIR/tmp/snv_intersect_unsorted.vcf.gz"
tabix -p vcf $OUTDIR/snv.vcf.gz

#########################
echo "PREP FOR SV"
#########################

## Remove SVTRYPE==BND
bcftools view -e 'INFO/SVTYPE="BND"' -Oz -o $OUTDIR/sv.vcf.gz $SV_DIR/filter/sv_filter2.vcf.gz
tabix -p vcf $OUTDIR/sv.vcf.gz

## Convert sample name in VCF
bcftools query -l $OUTDIR/sv.vcf.gz | sort > $OUTDIR/sv_samples.txt
bcftools reheader -s $RENAME $OUTDIR/sv.vcf.gz -o $OUTDIR/sv_renamed.vcf.gz

###############################
echo "Combining SV+SNV"
###############################
bcftools query -l $OUTDIR/sv_renamed.vcf.gz | sort > $OUTDIR/sv_samples.txt
bcftools query -l $OUTDIR/snv.vcf.gz | sort > $OUTDIR/snv_samples.txt
comm -12 $OUTDIR/snv_samples.txt $OUTDIR/sv_samples.txt > $OUTDIR/common_samples.txt

## Subset VCFs to common samples
bcftools view -S $OUTDIR/common_samples.txt --min-ac 1 --threads 4 -Oz -o $OUTDIR/sv_common_sample.vcf.gz $OUTDIR/sv_renamed.vcf.gz
bcftools view -S $OUTDIR/common_samples.txt --min-ac 1 --threads 4 -Oz -o $OUTDIR/snv_common_sample_before_resetId.vcf.gz $OUTDIR/snv.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' $OUTDIR/snv_common_sample_before_resetId.vcf.gz -Oz -o $OUTDIR/snv_common_sample.vcf.gz

tabix -p vcf $OUTDIR/sv_common_sample.vcf.gz
tabix -p vcf $OUTDIR/snv_common_sample.vcf.gz

# Concat them together (bcftools concat requires sample names are in the same order. bcftools view -S common_samples.txt achieves this).
bcftools concat --allow-overlaps --rm-dups none -Ou $OUTDIR/snv_common_sample.vcf.gz $OUTDIR/sv_common_sample.vcf.gz --threads 4 | \
bcftools sort -Oz --temp-dir "$OUTDIR/tmp/" -o $OUTDIR/merged.vcf.gz
tabix -p vcf $OUTDIR/merged.vcf.gz

#######################################
echo "Convert VCF to PLINK binary format"
#######################################
# Convert VCF to PLINK binary format
plink2 --vcf $OUTDIR/merged.vcf.gz \
       --max-alleles 2 \
       --mac 2 \
       --make-pgen \
       --vcf-half-call m \
       --out $OUTDIR/merged

echo "Generate Kinship Coefficients"
plink2 --pfile $OUTDIR/merged --make-king bin triangle --out $OUTDIR/cohort_kinship
plink2 --pfile $OUTDIR/merged --king-cutoff $OUTDIR/cohort_kinship 0.0442 --out $OUTDIR/filtered_samples  # 0.0442–0.0884 <- 3rd degree rel.

echo "Calculate r-squared for a specific SV "hit" against all other variants"
bcftools view -H "$OUTDIR/sv_common_sample.vcf.gz" | awk '{print $3}' > "$OUTDIR/sv_id.txt"

plink2 --pfile $OUTDIR/merged \
       --keep $OUTDIR/filtered_samples.king.cutoff.in.id \
       --mac 3 \
       --r2-unphased \
       --ld-snp-list $OUTDIR/sv_id.txt \
       --ld-window-kb 1000 \
       --ld-window-r2 0.5 \
       --out $OUTDIR/results_ld