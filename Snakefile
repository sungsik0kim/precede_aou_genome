import pandas as pd

MERGED_SNV_PATH='/home/jupyter/workspaces/longreadseqcontrolset/output/precede_aou_tier1/snv/merged'
SV_JOINT_CALL_PATH='/home/jupyter/workspaces/longreadseqcontrolset/output/precede_aou_tier1/sniffles2/joint_call'

df_aou_manifest = pd.read_csv(
                config['ref']['lrWGS_manifest'],
                dtype=str
            )

tmpdir= config['tmpdir'] 
outdir= config['outdir']
aou_snf_dir= config['input_data']['aou_snf_dir']
cadd_sv_dir = config['cadd_sv_dir']

ref_genome = config['ref']['genome_build']
ref_genome_fa = config['ref']['genome_fasta']
tandem_repeat = config['ref']['tandem_repeats_hg38']

snv_vcf = {
    'precede' : config['snv_vcf']['precede_filtered']
    ,'v7' : config['snv_vcf']['v7']
    ,'v8_BCM_Rev_high' : config['snv_vcf']['v8_BCM_Rev_high']
    ,'v8_BCM_Seq_high' : config['snv_vcf']['v8_BCM_Seq_high']
    ,'v8_BI_Rev_mid' : config['snv_vcf']['v8_BI_Rev_mid']
    ,'v8_BI_Seq_high' : config['snv_vcf']['v8_BI_Seq_high']
    ,'v8_BI_Seq_mid' : config['snv_vcf']['v8_BI_Seq_mid']
    ,'v8_HA_Rev_mid' : config['snv_vcf']['v8_HA_Rev_mid']
    ,'v8_UW_Rev_high' : config['snv_vcf']['v8_UW_Rev_high']
    ,'v8_UW_Seq_high' : config['snv_vcf']['v8_UW_Seq_high']
}
        
# utility funcs
def _get_samples():
    return df_aou_manifest['research_id'].unique()

def _get_file(tbl_column, index=False):
    ext = '.bai' if 'bam' in tbl_column else ''
    def f(wcs):
        value = df_aou_manifest\
                    .query('research_id == @wcs.sample')\
                    [tbl_column].iloc[0]
        if index:
            return f'{value}{ext}'
        else:
            return value
    return f

target_sv = []
target_snv = []
target_dev = []

target_sv.extend(
    [
        f'{outdir}/result/sv_prioritize.csv'
    ]
)

CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']
target_snv.extend(
    [
        f'{outdir}/result/snv_prioritize.csv'
    ]
)

target_dev.extend(
    [
        
    ]
)

bam_cols={
        "hg38":"grch38_bam"
        ,"t2t":"chm13v2.0_bam"
    }
aou_manifest_bam = bam_cols['hg38']
shell.prefix("set +u; ")

rule all:
    input: 
        f'{outdir}/result/snv_prioritize.csv',
        f'{outdir}/result/sv_prioritize.csv',
        f'{outdir}/result/phewas_variants.txt'
        

rule sv:
    input: 
        target_sv
        
rule snv:
    input: 
        target_snv

rule dev:
    input:
        target_dev
        
rule copy_bam:
    output:
        bam= temp(f'{tmpdir}/{{sample}}.{ref_genome}.bam')
    params:
        source= _get_file(aou_manifest_bam),
        snf= f'{outdir}/sniffles2/call/{{sample}}.sniffles2.{ref_genome}.snf'
    log:
        f'{outdir}/sniffles2/log/{{sample}}.log'
    threads: 1
    resources:
        copyq=1,
        mem_mb=2000
    shell:
        '''
        gsutil -u $GOOGLE_PROJECT cp {params.source} {output.bam}
        '''

rule index_bam:
    input:
        bam= f'{tmpdir}/{{sample}}.{ref_genome}.bam'
    output:
        bai= f'{tmpdir}/{{sample}}.{ref_genome}.bam.bai'
    threads: 4
    resources:
         mem_mb=2000
    priority: 50
    group: "sniffles_chain"
    shell:
        '''
        samtools index -@ {threads} {input.bam}
        '''     

####################################################################################
####### SV #########################################################################
####################################################################################

rule sniffles_call:
    input:
        bam= f'{tmpdir}/{{sample}}.{ref_genome}.bam',
        bai= f'{tmpdir}/{{sample}}.{ref_genome}.bam.bai'
    output:
        snf= f'{aou_snf_dir}/{{sample}}.sniffles2.{ref_genome}.snf',
        vcf= f'{aou_snf_dir}/{{sample}}.sniffles2.{ref_genome}.vcf.gz',
        idx= f'{aou_snf_dir}/{{sample}}.sniffles2.{ref_genome}.vcf.gz.tbi'
    params:
        ref= ref_genome_fa,
        trf= tandem_repeat
    log:
        f'{outdir}/sniffles2/log/{{sample}}.log'
    threads: 8
    resources:
         mem_mb=8000
    priority: 100
    group: "sniffles_chain"
    shell:
        '''
        sniffles \
            -i {input.bam} \
            --threads {threads} \
            --reference {params.ref} \
            --tandem-repeats {params.trf} \
            --phase \
            --minsvlen 50 \
            --snf {output.snf} \
            -v {output.vcf} \
            > {log} 2>&1
        '''
            
rule build_cohort:
    input:
        aou_snf_dir= config['input_data']['aou_snf_dir'],
        precede_snf_dir= config['input_data']['precede_snf_dir'],
        precede_mastersheet = config['ref']['precede_masterspreadsheet']
    output:
        ctrl_sample= f'{outdir}/cohort/ctrl_sample.csv',
        case_sample= f'{outdir}/cohort/case_sample.csv',
        internal_ctrl_sample= f'{outdir}/cohort/internal_ctrl_sample.csv',
        unaffected_rel_sample= f'{outdir}/cohort/unaffected_rel_sample.csv',
    log:
        f'{outdir}/cohort/log.txt'
    threads:
        2
    resources:
         mem_mb=4000
    params:
        script = config["script"]["build_cohort"],
        outdir = f'{outdir}/cohort'
    shell:
        """
        python {params.script} --aou_snf {input.aou_snf_dir} --precede_snf {input.precede_snf_dir} --precede_master {input.precede_mastersheet} --outdir {params.outdir} > {log} 2>&1
        """
            
rule sniffles_joint_call:
    input:
        ctrl_sample= f'{outdir}/cohort/ctrl_sample.csv',
        case_sample= f'{outdir}/cohort/case_sample.csv',
        internal_ctrl_sample= f'{outdir}/cohort/internal_ctrl_sample.csv',
        unaffected_rel_sample= f'{outdir}/cohort/unaffected_rel_sample.csv',
    log:
        f"{outdir}/sniffles2/joint_call/joint_call_log.txt"
    threads:
        24
    resources:
         mem_mb=400000
    output:
        joint_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.vcf.gz',
    params:
        vcf_unsorted= f"{tmpdir}/joint.unsorted.vcf",
        snf_list= f"{tmpdir}/snf_list.tsv",
        tmp= f'{outdir}/tmp'
    shell:
        """
        source /home/jupyter/.bashrc
        conda activate sniffles_2_6_3
        
        tail -n +2 {input.ctrl_sample} | cut -d',' -f2 > {params.snf_list}
        tail -n +2 {input.case_sample} | cut -d',' -f2 >> {params.snf_list}
        tail -n +2 {input.internal_ctrl_sample} | cut -d',' -f2 >> {params.snf_list}
        tail -n +2 {input.unaffected_rel_sample} | cut -d',' -f2 >> {params.snf_list}
        
        sniffles --input {params.snf_list} --vcf {params.vcf_unsorted} --threads {threads} --no-sort --allow-overwrite \
        > >(tee -a {log}) 2> >(tee -a {log} >&2)
        bcftools sort {params.vcf_unsorted} -Oz -o {output.joint_vcf} --temp-dir {params.tmp}  >> {log} 2>&1
        tabix -p vcf {output.joint_vcf}  >> {log} 2>&1
        """

rule make_target_bed:
    input:
        target_gene= config['input_data']['target_gene'],
        ld_block = config['input_data']['target_ld_block'],
        
    log:
        f"{outdir}/target_info/log.txt"
    threads:
        2
    resources:
        mem_mb=4000
    output:
        f'{outdir}/target_info/target_gene.sorted.bed.gz',
        f'{outdir}/target_info/target_gene_promoter.sorted.bed',
        f'{outdir}/target_info/target_gene_promoter.sorted.bed.gz',
        f'{outdir}/target_info/target_gene_target_ld.bed',
        f'{outdir}/target_info/target_gene_target_ld.bed.gz'
    params:
        outdir= f"{outdir}/target_info",
        expand_gene_list = config["script"]["expand_gene_list"],
        target_info= config["script"]["target_info"],
        target_gene_expanded = f'{outdir}/target_info/target_gene_expanded.txt',
        padding = 200000
    shell:
        """
        set -e -o pipefail
        cd {params.outdir}
        wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
        zcat gencode.v45.annotation.gtf.gz | grep -w "gene_name" > gencode.v45.annotation.gtf
        grep -w "gene" gencode.v45.annotation.gtf  > gene.gtf
        grep -w "transcript" gencode.v45.annotation.gtf | grep 'MANE_Select' > transcript.gtf
        grep -w "exon" gencode.v45.annotation.gtf | grep 'MANE_Select' > exon.gtf
        grep -w "UTR" gencode.v45.annotation.gtf | grep 'MANE_Select' > utr.gtf
        
        ## Append Gene list nearby LD and use this appended target gene list
        python {params.expand_gene_list} --target_ld {input.ld_block} --target_gene {input.target_gene} --transcript {params.outdir}/transcript.gtf --output {params.target_gene_expanded}
        
        python {params.target_info} --target {params.target_gene_expanded} --transcript {params.outdir}/transcript.gtf --padding {params.padding} --outdir {params.outdir}
        
        sort -k1,1V -k2,2n {params.outdir}/target_gene.bed > {params.outdir}/target_gene.sorted.bed
        bgzip -c {params.outdir}/target_gene.sorted.bed > {params.outdir}/target_gene.sorted.bed.gz
        tabix -p bed {params.outdir}/target_gene.sorted.bed.gz
        
        sort -k1,1V -k2,2n {params.outdir}/target_gene_promoter.bed > {params.outdir}/target_gene_promoter.sorted.bed
        bgzip -c {params.outdir}/target_gene_promoter.sorted.bed > {params.outdir}/target_gene_promoter.sorted.bed.gz
        tabix -p bed {params.outdir}/target_gene_promoter.sorted.bed.gz
        
        ## Combine with LD blocks
        cat <(tail -n +2 {input.ld_block} | awk 'BEGIN{{OFS="\t"}} {{print "chr"$1, $2, $3}}') \
            <(awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $3}}' {params.outdir}/target_gene.sorted.bed) \
            | sort -k1,1V -k2,2n | bedtools merge -i - > {params.outdir}/target_gene_target_ld.bed
        
        ## 2. Append the manual entry (non-MANE transcript)
        echo -e "chr3\t169564610\t169965060" >> {params.outdir}/target_gene_target_ld.bed
        sort -k1,1V -k2,2n {params.outdir}/target_gene_target_ld.bed -o {params.outdir}/target_gene_target_ld.bed
        
        bgzip -c {params.outdir}/target_gene_target_ld.bed > {params.outdir}/target_gene_target_ld.bed.gz
        tabix -p bed {params.outdir}/target_gene_target_ld.bed.gz
        """
        
rule intersect_vcf:
    input:
        joint_vcf= f'{SV_JOINT_CALL_PATH}/sniffles2_joint_call.{ref_genome}.vcf.gz',
        bed= f"{outdir}/target_info/target_gene_target_ld.bed.gz",
    log:
        f"{outdir}/sniffles2/filter_vcf_log.txt"
    threads:
        1
    resources:
        mem_mb=4000
    output:
        joint_vcf_filt= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.vcf.gz'
    params:
        outdir = f'{outdir}/sniffles2/joint_call'
    shell:
        """
        # Intersect: keep records that overlap BED intervals
        bcftools view -R {input.bed} -Oz -o {output.joint_vcf_filt} {input.joint_vcf}
        tabix -p vcf {output.joint_vcf_filt}
        """

rule run_truvari:
    input:
        joint_vcf_filt= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.vcf.gz',
        ctrl_sample= f'{outdir}/cohort/ctrl_sample.csv',
        case_sample= f'{outdir}/cohort/case_sample.csv',
        internal_ctrl_sample= f'{outdir}/cohort/internal_ctrl_sample.csv',
        unaffected_rel_sample= f'{outdir}/cohort/unaffected_rel_sample.csv',
    log:
        f"{outdir}/sniffles2/truvari_log.txt"
    threads:
        4
    resources:
        mem_mb=64000
    output:
        truvari_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari.vcf.gz',
        internal_ctrl_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari_internal_ctrl.vcf.gz',
        unaffected_rel_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari_unaffected_rel.vcf.gz',
    params:
        outdir = f'{outdir}/sniffles2/joint_call',
        case_ctrl_ids= f'{outdir}/sniffles2/joint_call/case_ctrl_ids.txt',
        internal_ctrl_ids= f'{outdir}/sniffles2/joint_call/internal_ctrl_ids.txt',
        unaffected_rel_ids= f'{outdir}/sniffles2/joint_call/unaffected_rel_ids.txt',
        ref = ref_genome_fa,
        tmp = f'{outdir}/tmp'
        
    shell:
        """
        source /home/jupyter/.bashrc
        conda activate truvari
        truvari collapse -i {input.joint_vcf_filt} \
            -o {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf \
            -c {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_collapsed.vcf \
            -f {params.ref} -k common
            
        bcftools sort {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf \
            -Oz -o {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf.gz \
            --temp-dir {params.tmp} 
        tabix -p vcf {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf.gz
        
        ## Split truvari output (1) AoU and Precede FPC (2) Internal control (3) Unaffected relatives 
        tail -n +2 {input.ctrl_sample} | cut -d',' -f1 | sed 's/$/.sniffles2.hg38/' > {params.case_ctrl_ids}
        tail -n +2 {input.case_sample} | cut -d',' -f1 >> {params.case_ctrl_ids}
        tail -n +2 {input.internal_ctrl_sample} | cut -d',' -f1 > {params.internal_ctrl_ids}
        tail -n +2 {input.unaffected_rel_sample} | cut -d',' -f1 > {params.unaffected_rel_ids}
        
        bcftools view -S {params.internal_ctrl_ids} --threads {threads} \
            {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf.gz \
            -Oz -o {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_internal_ctrl.vcf.gz &
        bcftools view -S {params.unaffected_rel_ids} --threads {threads} \
            {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf.gz \
            -Oz -o {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_unaffected_rel.vcf.gz &
        bcftools view -S {params.case_ctrl_ids} --threads {threads} \
            {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_beforeSplit.vcf.gz \
            -Oz -o {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari.vcf.gz &
        wait
        
        tabix -p vcf {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_internal_ctrl.vcf.gz &
        tabix -p vcf {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari_unaffected_rel.vcf.gz &
        tabix -p vcf {params.outdir}/sniffles2_joint_call.{ref_genome}.target_region.truvari.vcf.gz &
        wait
        """

rule sv_filter_1:
    input:
        vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari.vcf.gz',
        ctrl_sample= f'{outdir}/cohort/ctrl_sample.csv',
        case_sample= f'{outdir}/cohort/case_sample.csv'
    
    output:
        vcf= f'{outdir}/sv/filter/sv_filter1.vcf.gz',
        callset= f'{outdir}/sv/filter/sv_filter1.callset_only.vcf.gz',
    threads:
        4
    resources:
        mem_mb=48000
    params:
        script = config["script"]["sv_filter_1"],
        outdir = f'{outdir}/sv/filter',
        
        
    conda:
        "envs/analysis.yaml"
    log:
        f"{outdir}/sv/filter/sv_filter_1.log"   
    shell:
        """
        mkdir -p {params.outdir}
        
        ###
        tail -n +2 {input.case_sample} | cut -d',' -f1 > {params.outdir}/case_sample.txt
        tail -n +2 {input.ctrl_sample} | cut -d',' -f1 | sed 's/$/.sniffles2.hg38/' > {params.outdir}/case_ctrl_sample.txt
        cat {params.outdir}/case_sample.txt >> {params.outdir}/case_ctrl_sample.txt 

        # Filter VCF with ALT in case and get ID
        bcftools view -S {params.outdir}/case_sample.txt --force-samples --min-ac 1 --threads 4 -Oz -o {params.outdir}/case_only.vcf.gz {input.vcf} > {log}
        tabix -p vcf {params.outdir}/case_only.vcf.gz >>{log}
        bcftools query -f '%ID\n' {params.outdir}/case_only.vcf.gz > {params.outdir}/filt_variant.id

        # Filter the original vcf to get the final 
        bcftools view -S {params.outdir}/case_ctrl_sample.txt --force-samples -i 'ID=@{params.outdir}/filt_variant.id' -f PASS --threads 4 -Oz -o {output.vcf} {input.vcf} >>{log}
        tabix -p vcf {output.vcf}
        
        ###
        
        # Filtered VCF (Only callset for downstream operation)
        bcftools view --drop-genotypes {output.vcf} -Oz -o {output.callset}
        tabix -p vcf {output.callset}
        """
        
        
rule annotate_gnomAD:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter1.callset_only.vcf.gz',
        bed= config['ref']['svafotate_bed']
    output:
        f'{outdir}/sv/annotate/sv_filter1.gnomadAF.vcf'
    threads:
        4
    resources:
         mem_mb=48000
    conda:
        "envs/svafotate.yaml"
    log:
        f"{outdir}/sniffles2/joint_call/svafotate.log"   
    shell:
        """
        svafotate annotate -v {input.vcf} -o {output} -b {input.bed} -f 0.5 -a pops -s gnomAD --cpu {threads} > {log} 2>&1
        ## CHANGE : -a EUR -s gnomAD
        """

rule annotate_CoLoRSdb:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter1.callset_only.vcf.gz',
        CoLoRdb = config['ref']['CoLoRdb_sv']
    output:
        vcf= f'{outdir}/sv/annotate/sv_filter1.colorsAF.vcf'
    threads:
        2
    resources:
        mem_mb=24000
    shell:
        """
        svdb --query \
          --query_vcf  {input.vcf} \
          --db {input.CoLoRdb} \
          --overlap 0.5 \
          --bnd_distance 500 \
          --in_occ AC \
          --in_frq AF \
          --out_occ CoLoRdb_AC \
          --out_frq CoLoRdb_AF \
          > {output.vcf}
        """
        
rule annotate_hprc:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter1.callset_only.vcf.gz',
        hprc = config['ref']['hprc_sv']
    output:
        vcf= f'{outdir}/sv/annotate/sv_filter1.hprcAF.vcf'
    threads:
        2
    resources:
        mem_mb=24000
    shell:
        """
        svdb --query \
          --query_vcf  {input.vcf} \
          --db {input.hprc} \
          --overlap 0.5 \
          --bnd_distance 500 \
          --in_occ AC \
          --in_frq AF \
          --out_occ hprc_AC \
          --out_frq hprc_AF \
          > {output.vcf}
        """
        
rule merge_AF:
    input:
        gnomAD= f'{outdir}/sv/annotate/sv_filter1.gnomadAF.vcf',
        CoLoRdb = f'{outdir}/sv/annotate/sv_filter1.colorsAF.vcf',
        hprc = f'{outdir}/sv/annotate/sv_filter1.hprcAF.vcf'
    output:
        f'{outdir}/sv/annotate/sv_filter1.AF.csv'
    threads:
        2
    resources:
        mem_mb=24000
    params:
        script = config["script"]["merge_AF"],
    shell:
        """
        python {params.script} --gnomAD {input.gnomAD} --CoLoRdb {input.CoLoRdb} --hprc {input.hprc} --output {output}
        """

rule sv_filter_2:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter1.vcf.gz',
        af = f'{outdir}/sv/annotate/sv_filter1.AF.csv',
        case= f'{outdir}/cohort/case_sample.csv',
        ctrl= f'{outdir}/cohort/ctrl_sample.csv',
        internal_ctrl_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari_internal_ctrl.vcf.gz',
        unaffected_rel_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari_unaffected_rel.vcf.gz',
        master= config['ref']['precede_masterspreadsheet'],
        aou_covar = config['ref']['aou_covariate'],
        
    output:
        vcf= f'{outdir}/sv/filter/sv_filter2.vcf.gz',
        filt_csv = f'{outdir}/sv/filter/sv_filter2.csv',
        callset= f'{outdir}/sv/filter/sv_filter2.callset_only.vcf.gz'
        
    threads:
        16
    resources:
        mem_mb=48000
    params:
        script = config["script"]["sv_filter_2"],
        outdir = f'{outdir}/sv/filter',
        filt_bed = f'{outdir}/sv/filter/sv_filter2.bed',
        id= f'{outdir}/sv/filter/sv_filter2.id',
        
    conda:
        "envs/analysis.yaml"
    log:
        f"{outdir}/sv/filter/sv_filter_2.log"   
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} --vcf {input.vcf} --af {input.af} --case {input.case} --ctrl {input.ctrl} \
            --master {input.master} --aou_covar {input.aou_covar} --ctrl2 {input.internal_ctrl_vcf} --ctrl3 {input.unaffected_rel_vcf} \
            --outdir {params.outdir} > {log} 2>&1
        
        # Filtered VCF
        bcftools view -R {params.filt_bed} -i 'ID=@{params.id}' --threads {threads} {input.vcf} -Oz  -o {output.vcf}
        tabix -p vcf {output.vcf}
        
        # Filtered VCF (Only callset for downstream operation)
        bcftools view --drop-genotypes {output.vcf} -Oz -o {output.callset}
        tabix -p vcf {output.callset}
        """

rule find_surrogate_snp: # Some code blocks are commented out. BE CAUTIOUS
    input:
        sv= f'{outdir}/sv/filter/sv_filter2.vcf.gz',
        snv= expand(f'{MERGED_SNV_PATH}/AoU_Precede_chr{{chrom}}.vcf.gz', chrom=CHROMS),
        sniffles_rename = config["ref"]["sniffles_rename"]
            
    output:
        f'{outdir}/plink/surrogate_snp.tsv'
    params:
        sv_dir = f'{outdir}/sv',
        snv_dir = f'{MERGED_SNV_PATH}',
        out_dir = f'{outdir}/plink',
        find_surrogate_snp_sh = config["script"]["find_surrogate_snp_sh"],
        find_surrogate_snp_py = config["script"]["find_surrogate_snp_py"],
    threads:
        48
    resources:
        mem_mb=256000
    log: 
        f'{outdir}/log/find_surrogate_snp.log'
    shell:
        """
        bash {params.find_surrogate_snp_sh} {params.sv_dir} {params.snv_dir} {params.out_dir} {input.sniffles_rename} > {log} 2>&1
        python {params.find_surrogate_snp_py} --plink_ld {params.out_dir}/results_ld.vcor --outcsv {output} >> {log} 2>&1
        """        
        
rule vep:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter2.callset_only.vcf.gz'
    log:
        f"{outdir}/sv/annotate/vep.log"
    threads:
        8
    output:
        vep = f'{outdir}/sv/annotate/sv_filter2.vep.vcf'
    params:
        script = config["script"]["vep_sv"],
        outdir = f'{outdir}/sv/annotate',
        prefix = 'sv_filter2'
    shell:
        """
        sh {params.script} {input.vcf} {params.outdir} {params.prefix} > {log} 2>&1
        """


rule prep_cadd_sv:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter2.callset_only.vcf.gz'
    output:
        bed_temp= temp(f'{cadd_sv_dir}/input/temp_id_precede.bed'),
        bed = temp(f'{cadd_sv_dir}/input/id_precede.bed')
    threads:
        2
    resources:
        mem_mb=8000
    params:
        script = config["script"]["prep_cadd_sv"] 
    log: 
        f'{outdir}/sv/annotate/prep_cadd_sv.log'
    
    shell:
        """
        python {params.script} --vcf_in {input.vcf} --bed_out {output.bed_temp} > {log} 2>&1
        bedtools sort -i {output.bed_temp} > {output.bed}
        """
        
rule cadd_sv:
    input:
        bed= f'{cadd_sv_dir}/input/id_precede.bed'
    output:
        bed_temp= temp(f'{cadd_sv_dir}/output/precede_score.bed'),
        bed= f'{outdir}/sv/annotate/sv_filter2.CADD.bed'
    threads:
        8
    resources:
        mem_mb=48000
    log: 
        f'{outdir}/sv/annotate/cadd_sv.log'
    shell:
        """
        source /home/jupyter/.bashrc
        conda activate run.caddsv
        
        cd {cadd_sv_dir}
        echo 'dataset: precede' > config.yml 
        snakemake --use-conda --rerun-incomplete --configfile config.yml -j {threads} > {log} 2>&1
        cp {output.bed_temp} {output.bed}
        """
        
rule AnnotSV:
    input:
        vcf= f'{outdir}/sv/filter/sv_filter2.vcf.gz',
    output:
        tsv= f'{outdir}/sv/annotate/sv_filter2.AnnotSV.tsv',
    params:
        annotsv= config["AnnotSV_dir"]
    threads:
        2
    resources:
        mem_mb=16000
    log: 
        f'{outdir}/sv/annotate/annotSV.log'
    shell:
        """
        {params.annotsv}/bin/AnnotSV \
          -SVinputFile {input.vcf} \
          -outputFile {output.tsv} \
          -genomeBuild GRCh38 \
          -overlap 50 \
          -reciprocal 1 \
          -promoterSize 1000 \
          -full_annotation > {log} 2>&1
        """
        
rule family_segregation :
    input:
        vcf= f'{outdir}/sv/filter/sv_filter2.vcf.gz',
        unaffected_rel_vcf= f'{outdir}/sniffles2/joint_call/sniffles2_joint_call.{ref_genome}.target_region.truvari_unaffected_rel.vcf.gz',
        master= config['ref']['precede_masterspreadsheet'],
    output:
        f'{outdir}/sv/annotate/sv_filter2.family_segregation.csv'
    params:
        script= config["script"]["family_segregation"]
    threads:
        2
    resources:
        mem_mb=16000
    log: 
        f'{outdir}/sv/annotate/family_segregation.log'
    shell:
        """
        python {params.script} --vcf {input.vcf} --unaffected_rel_vcf {input.unaffected_rel_vcf} --master {input.master} --output {output}
        """

rule sv_prioritize:
    input:
        filt2_csv= f'{outdir}/sv/filter/sv_filter2.csv',
        af= f'{outdir}/sv/annotate/sv_filter1.AF.csv',
        vep= f'{outdir}/sv/annotate/sv_filter2.vep.vcf',
        cadd= f'{outdir}/sv/annotate/sv_filter2.CADD.bed',
        annotSV= f'{outdir}/sv/annotate/sv_filter2.AnnotSV.tsv',
        family_segregation = f'{outdir}/sv/annotate/sv_filter2.family_segregation.csv',
        master= config['ref']['precede_masterspreadsheet']
            
    output:
        f'{outdir}/result/sv_prioritize.csv'
    params:
        script = config["script"]["sv_prioritize"]
    threads:
        2
    resources:
        mem_mb=8000
    log: 
        f'{outdir}/log/sv_prioritize.log'
    shell:
        """
        python {params.script} \
            --filt2_csv {input.filt2_csv} \
            --af {input.af} \
            --vep {input.vep} \
            --cadd {input.cadd} \
            --annotSV {input.annotSV} \
            --family_segregation {input.family_segregation} \
            --master {input.master} \
            --outcsv {output} \
            > {log} 2>&1
        """

####################################################################################
####### SNV ########################################################################
####################################################################################
rule qualFilter_precede_snv:
    input:
        config['snv_vcf']['precede_unfiltered']
    output:
        config['snv_vcf']['precede_filtered']
    threads:
        8 
    params:
        ref= ref_genome_fa
    shell:
        """
        bcftools filter \
            -e 'QUAL < 40' \
            -s 'lowQual' \
            -m + \
            -O z \
            --threads {threads} \
            -o {output} \
            {input}
        """

rule normalize_snv:
    input:
        vcf= lambda wildcards: snv_vcf[wildcards.batch]
    output:
        vcf_norm=f'{outdir}/snv/normalized/{{batch}}.vcf.unsorted.gz'
    threads:
        8 
    params:
        ref= ref_genome_fa
    shell:
        """
        bcftools norm --threads {threads} -m -any -f {params.ref} {input.vcf} -O z -o {output.vcf_norm}
        """

rule sort_index_snv:
    input:
        vcf_norm= f'{outdir}/snv/normalized/{{batch}}.vcf.unsorted.gz'
    output:
        vcf_norm_sorted=f'{outdir}/snv/normalized/{{batch}}.vcf.gz'
    threads:
        1 
    params:
        ref= ref_genome_fa,
        tmp= f'{tmpdir}/{{batch}}'
    shell:
        """
        bcftools sort {input.vcf_norm} -O z -o {output.vcf_norm_sorted} --temp-dir {params.tmp}
        tabix -p vcf {output.vcf_norm_sorted}
        rm {input.vcf_norm}
        """

rule merge_snv:
    input:
        vcfs = expand(f'{outdir}/snv/normalized/{{batch}}.vcf.gz', batch=snv_vcf.keys())
    output:
        vcf_merged=f'{outdir}/snv/merged/AoU_Precede_chr{{chrom}}.vcf.gz'
    log: f'{outdir}/snv/merged/log_chr{{chrom}}.txt'
    params:
        region= f'chr{{chrom}}'
    threads:
        8
    shell:
        """
        bcftools merge \
            --missing-to-ref \
            --force-samples \
            --threads {threads} \
            -m none \
            -r {params.region} \
            -O z \
            -o {output.vcf_merged} \
            {input.vcfs} > {log} 2>&1
        tabix -p vcf {output.vcf_merged}
        """
        
        
rule intersect_vcf_snv:
    input:
#         vcf= f'{outdir}/snv/merged/AoU_Precede_chr{{chrom}}.vcf.gz',
        vcf= f'{MERGED_SNV_PATH}/AoU_Precede_chr{{chrom}}.vcf.gz',
        bed = f'{outdir}/target_info/target_gene_target_ld.bed.gz'
    log:
        f"{outdir}/snv/target/intersect_chr{{chrom}}_vcf_snv_log.txt"
    threads:
        2
    resources:
        mem_mb=8000
    output:
        vcf= f'{outdir}/snv/target/AoU_Precede_chr{{chrom}}_snv_target.vcf.gz'
    params:
        outdir = f'{outdir}/snv/target'
    shell:
        """
        # Intersect: keep records that overlap BED intervals
        bcftools view --threads {threads} -R {input.bed} -Oz -o {output.vcf} {input.vcf}
        tabix -p vcf {output.vcf}
        """
        
rule combine_vcf_chr:
    input:
        vcfs= expand(f'{outdir}/snv/target/AoU_Precede_chr{{chrom}}_snv_target.vcf.gz', chrom=CHROMS),
    log:
        f'{outdir}/snv/target/combine_vcf_chr_log.txt'
    threads:
        4
    resources:
        mem_mb=8000
    output:
        vcf= f'{outdir}/snv/target/AoU_Precede_snv_target.vcf.gz',
        temp_vcf = temp(f'{outdir}/snv/target/AoU_Precede_snv_target_temp.vcf.gz')
    params:
        outdir = f'{outdir}/snv/target',
        tmp= f'{outdir}/tmp'
    shell:
        """
        bcftools concat {input.vcfs} -Oz -o {output.temp_vcf}
        bcftools sort {output.temp_vcf} --temp-dir {params.tmp} -Oz -o {output.vcf} 
        tabix -p vcf {output.vcf}
        """

rule snv_filter_1:
    input:
        vcf= f'{outdir}/snv/target/AoU_Precede_snv_target.vcf.gz',
        ctrl_sample= f'{outdir}/cohort/ctrl_sample.csv',
        case_sample= f'{outdir}/cohort/case_sample.csv'
    output:
        vcf = f'{outdir}/snv/filter/snv_filter1.vcf.gz',
        id = f'{outdir}/snv/filter/filt_variant.id'
    log:
        f'{outdir}/snv/filter/snv_filter_1.log'
    threads:
        4
    resources:
        mem_mb=48000
    params:
        outdir = f'{outdir}/snv/filter',
        script = config["script"]["snv_filter_1"],
    shell:
        """
        mkdir -p {params.outdir}
        tail -n +2 {input.case_sample} | cut -d',' -f1 > {params.outdir}/case_sample.txt
        tail -n +2 {input.ctrl_sample} | cut -d',' -f1 > {params.outdir}/case_ctrl_sample.txt
        cat {params.outdir}/case_sample.txt >> {params.outdir}/case_ctrl_sample.txt 

        # Filter VCF with ALT in case and get ID
        bcftools view -S {params.outdir}/case_sample.txt --force-samples --min-ac 1 --threads 4 -Oz -o {params.outdir}/case_only.vcf.gz {input.vcf} > {log}
        tabix -p vcf {params.outdir}/case_only.vcf.gz >>{log}
        bcftools query -f '%ID\n' {params.outdir}/case_only.vcf.gz > {params.outdir}/_filt_variant.id

        # Remove >=3 alt alleles 
        grep -v ";.*;" {params.outdir}/_filt_variant.id > {output.id}

        # Filter the original vcf to get the final 
        bcftools view -S {params.outdir}/case_ctrl_sample.txt  --force-samples -i 'ID=@{output.id}' -f PASS --threads 4 -Oz -o {output.vcf} {input.vcf} >>{log}
        tabix -p vcf {output.vcf}
        
        # Clean-up
        rm {params.outdir}/case_only.vcf.gz >>{log}
        rm {params.outdir}/case_only.vcf.gz.tbi >>{log}
        rm {params.outdir}/_filt_variant.id >>{log}
        """
        
rule prep_precede_ctrl:
    input:
        vcf= f'{outdir}/snv/target/AoU_Precede_snv_target.vcf.gz',
        id= f'{outdir}/snv/filter/filt_variant.id'
    output:
        internal_ctrl_vcf_snv=f'{outdir}/snv/target/AoU_Precede_snv_target.internal_ctrl.vcf.gz',
        unaffected_rel_vcf_snv=f'{outdir}/snv/target/AoU_Precede_snv_target.unaffected_rel.vcf.gz',
    log:
        f'{outdir}/log/prep_precede_ctrl.log'
    threads:
        4
    resources:
        mem_mb=16000
    params:
        internal_ctrl_ids= f'{outdir}/sniffles2/joint_call/internal_ctrl_ids.txt',
        unaffected_rel_ids= f'{outdir}/sniffles2/joint_call/unaffected_rel_ids.txt',
    shell:
        """
        bcftools view -S {params.internal_ctrl_ids} -i 'ID=@{input.id}' --threads {threads} {input.vcf} \
            -Oz -o {output.internal_ctrl_vcf_snv} &
        bcftools view -S {params.unaffected_rel_ids} -i 'ID=@{input.id}' --threads {threads} {input.vcf} \
            -Oz -o {output.unaffected_rel_vcf_snv} &
        wait
        
        tabix -p vcf {output.internal_ctrl_vcf_snv} &
        tabix -p vcf {output.unaffected_rel_vcf_snv} &
        wait
        """

rule snv_annotate_af:
    input:
        vcf= f'{outdir}/snv/filter/snv_filter1.vcf.gz',
        hprc_vcf= config['ref']['hprc'],
        CoLoRdb_vcf= config['ref']['CoLoRdb_snv']
    output:
        af= f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt_anno_AF.vcf.gz',
    threads:
        16
    resources:
        mem_mb=4000
    params:
        outdir= f'{outdir}/snv/annotate'
    log: 
        f'{outdir}/snv/annotate/snv_annotate_af.log',
    shell:
        """
        bcftools view --drop-genotypes {input.vcf} -Oz -o {params.outdir}/temp.vcf.gz
        tabix -p vcf {params.outdir}/temp.vcf.gz
        
        bcftools annotate \
              --threads {threads} \
              -a {input.hprc_vcf} \
              -c INFO/HPRC_AF:=INFO/AF \
              -Oz -o {params.outdir}/temp2.vcf.gz {params.outdir}/temp.vcf.gz
        tabix -p vcf {params.outdir}/temp2.vcf.gz
              
        bcftools annotate \
              --threads {threads} \
              -a {input.CoLoRdb_vcf} \
              -c INFO/COLOR_AF:=INFO/AF \
              -Oz -o {output.af} {params.outdir}/temp2.vcf.gz
              
        tabix -p vcf {output.af}
        rm {params.outdir}/temp.vcf.gz {params.outdir}/temp.vcf.gz.tbi \
           {params.outdir}/temp2.vcf.gz {params.outdir}/temp2.vcf.gz.tbi
        """

rule snv_gnomad_af:
    input:
        vcf= f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt_anno_AF.vcf.gz', # To use genotype dropped vcf
    log:
        f"{outdir}/snv/annotate/snv_gnomad_af.log"
    threads:
        16
    resources:
        mem_mb=28000
    output:
        f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt.gnomad.vcf'
    params:
        script = config["script"]["vep_gnomad"],
        outdir = f'{outdir}/snv/annotate',
        prefix = f'AoU_Precede_snv_target_filt'
    shell:
        """
        sh {params.script} {input.vcf} {params.outdir} {params.prefix} > {log} 2>&1
        """
        
rule annotate_dbnsfp:
    input:
        vcf= f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt_anno_AF.vcf.gz', # To use genotype dropped vcf
    output:
        vcf = f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt.dbnsfp.vcf'
    conda:
        "envs/snpsift.yaml"
    log:
        f"{outdir}/snv/annotate/dbnsfp.log"
    threads:
        2
    resources:
        mem_mb=28000
    params:
        dbnsfp= config['ref']['dbnsfp']
    shell:
        """
        SnpSift dbnsfp -v -db {params.dbnsfp} -a -f Interpro_domain,REVEL_score,REVEL_rankscore,VEST4_score,VEST4_rankscore,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,ESM1b_score,ESM1b_converted_rankscore,ESM1b_pred,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,VEP_canonical,MANE,dbNSFP_POPMAX_AF {input.vcf} > {output.vcf} 2> {log}
        """
        
rule snv_filter_2:
    input:
        af= f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt_anno_AF.vcf.gz',
        gnomad = f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt.gnomad.vcf',
        dbnsfp = f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt.dbnsfp.vcf',
        vcf= f'{outdir}/snv/filter/snv_filter1.vcf.gz',
        case= f'{outdir}/cohort/case_sample.csv',
        ctrl= f'{outdir}/cohort/ctrl_sample.csv',
        internal_ctrl_vcf=f'{outdir}/snv/target/AoU_Precede_snv_target.internal_ctrl.vcf.gz',
        unaffected_rel_vcf=f'{outdir}/snv/target/AoU_Precede_snv_target.unaffected_rel.vcf.gz',
        master= config['ref']['precede_masterspreadsheet'],
        aou_covar = config['ref']['aou_covariate'],
    log:
        f"{outdir}/snv/filter/snv_filter_2.log"
    threads:
        8
    resources:
        mem_mb=28000
    output:
        csv= f'{outdir}/snv/annotate/snv_filter2.csv',
        id= f'{outdir}/snv/annotate/snv_filter2_id.txt',
        bed= f'{outdir}/snv/annotate/snv_filter2.bed',
        vcf= f'{outdir}/snv/annotate/snv_filter2.vcf.gz',
        callset= f'{outdir}/snv/annotate/snv_filter2.callset_only.vcf.gz',
    params:
        script = config["script"]["snv_filter_2"],
        outdir = f'{outdir}/snv/annotate'
    shell:
        """
        python {params.script} --vcf {input.vcf} \
                --af {input.af} \
                --gnomad {input.gnomad} \
                --dbnsfp {input.dbnsfp} \
                --outdir {params.outdir} \
                --case {input.case} \
                --ctrl {input.ctrl} \
                --ctrl2 {input.internal_ctrl_vcf} \
                --ctrl3 {input.unaffected_rel_vcf} \
                --master {input.master} \
                --aou_covar {input.aou_covar} \
                > {log} 2>&1
        bcftools view -R {output.bed} -i 'ID=@{output.id}' --threads {threads} {input.vcf} -Oz  -o {output.vcf}
        bcftools view --drop-genotypes {output.vcf} -Oz -o {output.callset}
        
        tabix -p {output.vcf} &
        tabix -p {output.callset} &
        wait
        
        """
        
rule family_segregation_snv :
    input:
        vcf= f'{outdir}/snv/annotate/snv_filter2.vcf.gz',
        unaffected_rel_vcf= f'{outdir}/snv/target/AoU_Precede_snv_target.unaffected_rel.vcf.gz',
        master= config['ref']['precede_masterspreadsheet'],
    output:
        f'{outdir}/snv/annotate/snv_filter2.family_segregation.csv'
    params:
        script= config["script"]["family_segregation"]
    threads:
        2
    resources:
        mem_mb=16000
    log: 
        f'{outdir}/snv/annotate/family_segregation.log'
    shell:
        """
        python {params.script} --vcf {input.vcf} --unaffected_rel_vcf {input.unaffected_rel_vcf} --master {input.master} --output {output}
        """

rule vep_snv:
    input:
        vcf= f'{outdir}/snv/annotate/snv_filter2.callset_only.vcf.gz', # To use genotype dropped vcf
    log:
        f"{outdir}/snv/annotate/vep_snv.log"
    threads:
        28
    resources:
        mem_mb=48000
    output:
        vep = f'{outdir}/snv/annotate/snv_filter2.vep.vcf'
    params:
        script = config["script"]["vep_snv"],
        outdir = f'{outdir}/snv/annotate',
        prefix = f'snv_filter2'
    shell:
        """
        sh {params.script} {input.vcf} {params.outdir} {params.prefix} > {log} 2>&1
        """

rule prep_promoterai:
    input:
        vcf= f'{outdir}/snv/annotate/snv_filter2.callset_only.vcf.gz', # To use genotype dropped vcf
        promoter_bed = f'{outdir}/target_info/target_gene_promoter.sorted.bed'
    output:
        vcf_norm = temp(f'{outdir}/snv/annotate/snv_filter2.norm.vcf.gz'),
        tsv = f'{outdir}/snv/annotate/promoterai_input.tsv'
    threads:
        8
    resources:
        mem_mb=26000
    params:
        script= config['script']['prep_promoterai'],
        ref= ref_genome_fa,
    log:
        tsv = f'{outdir}/snv/annotate/prep_promoterai.log'
    shell:
        """
        bcftools norm --threads {threads} -m -any -f {params.ref} {input.vcf} -O z -o {output.vcf_norm}
        python {params.script} --vcf {output.vcf_norm} --promoter_bed {input.promoter_bed} --output {output.tsv} > {log} 2>&1
        """

rule run_promoterai_GPU:
    input:
        tsv= f'{outdir}/snv/annotate/promoterai_input.tsv'
    output:
        tsv = f'{outdir}/snv/annotate/promoterai_input.promoterAI_v1_hg38_finetune.tsv'
    threads:
        8
    resources:
        mem_mb=26000
    params:
        model= config['promoterAI_model'],
        ref= ref_genome_fa,
    log:
        tsv = f'{outdir}/snv/annotate/run_promoterai.log'
    shell:
        """
        source /home/jupyter/.bashrc
        set +u
        conda activate tf-gpu
        
        promoterai \
            --model_folder {params.model} \
            --var_file {input.tsv} \
            --fasta_file {params.ref} \
            --input_length 20480  > {log} 2>&1
        """

rule snv_prioritize:
    input:
        filt2_csv= f'{outdir}/snv/annotate/snv_filter2.csv',
        dbnsfp= f'{outdir}/snv/annotate/AoU_Precede_snv_target_filt.dbnsfp.vcf',
        vep= f'{outdir}/snv/annotate/snv_filter2.vep.vcf',
        promoterai= f'{outdir}/snv/annotate/promoterai_input.promoterAI_v1_hg38_finetune.tsv',
        family_segregation = f'{outdir}/snv/annotate/snv_filter2.family_segregation.csv',
        master= config['ref']['precede_masterspreadsheet']
    output:
        f'{outdir}/result/snv_prioritize.csv'
    params:
        script = config["script"]["snv_prioritize"]
    threads:
        2
    resources:
        mem_mb=8000
    log: 
        f'{outdir}/log/snv_prioritize.log'
    shell:
        """
        python {params.script} \
            --filt2_csv {input.filt2_csv} \
            --dbnsfp {input.dbnsfp} \
            --vep {input.vep} \
            --promoterai {input.promoterai} \
            --family_segregation {input.family_segregation} \
            --master {input.master} \
            --outcsv {output} \
            > {log} 2>&1
        """
        
####################################################################################
####### Discover Putative Hits #####################################################
####################################################################################
rule prep_PheWAS:
    input:
        snv=f'{outdir}/result/snv_prioritize.csv',
        sv=f'{outdir}/plink/surrogate_snp.tsv'
    output:
        f'{outdir}/result/phewas_variants.txt'
    shell:
        """
        python3 -c "import sys, csv; reader = csv.DictReader(sys.stdin, delimiter='\\t'); [print(row['ID_B_2']) for row in reader]" < {input.sv} > {output}
        python3 -c 'import sys, pandas as pd; df = pd.read_csv(sys.stdin); sub = df[df["prioritized_variant"] == 1]; print("\\n".join(sub["CHROM"].str[3:] + "-" + sub["POS"].astype(str) + "-" + sub["REF"] + "-" + sub["ALT"]))' < {input.snv}  >> {output}
        """
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
