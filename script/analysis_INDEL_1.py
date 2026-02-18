import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.ticker import LogLocator
import numpy as np
import os
import argparse
from scipy.stats import chi2_contingency


def main(vep_csv, precede_mastersheet, outdir):
    ###################
    ## Read data
    ###################
    print("## Read data ##")
    df = pd.read_csv(vep_csv)
    df_precede_master = pd.read_csv(precede_mastersheet,dtype=str)
    
    ###################
    ## Subset samples
    ###################
    print("## Subset samples ##")
    # AoU
    aou_cols = [col for col in df.columns if 'sniffles2.hg38' in col] ################ [!!] FUTURE: Need to select cancer-history-free casees
    # Cohort
    mask = (df_precede_master['Exclusion']!='NOT CONFIDENT') & (df_precede_master['Affected Status']=='Affected')
    cohort_cols = df_precede_master.loc[mask, "Analysis_sample_id"].tolist()
    cohort_cols = list(set(cohort_cols) & set(df.columns))
    
    print("Usable AoU Cohort Size: {}".format(len(aou_cols)))
    print("Usable FPC Cohort size: {}".format(len(cohort_cols)))
    print("Total sample size in the master sheet: {}".format(df_precede_master.shape[0]))
    print("Affecte sample size in the master sheet: {}".format(sum(mask)))
    ###################
    ## Plot SV size distribution
    ###################
    print("## Plot SV size distribution ##")
    sv_min=50
    sv_max=20000
    interval=10
    
    bins = np.linspace(sv_min, sv_max, num=int(sv_max/interval))
    bins_neg = np.linspace(-sv_max, -sv_min, num=int(sv_max/interval))
    n_bin= len(bins)
    l_start = bins_neg[:n_bin-1].tolist()+ bins[:n_bin-1].tolist()
    l_end = bins_neg[1:n_bin].tolist()+ bins[1:n_bin].tolist()
    df_svlen_bin = pd.DataFrame({"start":l_start, "end":l_end})
    
    l_aou_count = []
    l_cohort_count = []
    
    for _,row in tqdm(df_svlen_bin.iterrows()):
        start, end = row['start'], row['end']
        mask = (df['SVLEN']>=start) & (df['SVLEN']<end)
        N=sum(mask)
        if N==0: 
            l_aou_count.append(np.nan)
            l_cohort_count.append(np.nan)
            continue
        else:
            aou_AC = df.loc[mask,aou_cols].sum().sum()
            aou_max = N*len(aou_cols)
            aou_na = df.loc[mask,aou_cols].isna().sum().sum()
            aou_AC_adj = aou_AC/((aou_max-aou_na)/aou_max)
            aou_AC_per_case = aou_AC_adj/len(aou_cols)
            l_aou_count.append(aou_AC_per_case)

            cohort_AC = df.loc[mask,cohort_cols].sum().sum()
            cohort_max = N*len(cohort_cols)
            cohort_na = df.loc[mask,cohort_cols].isna().sum().sum()
            cohort_AC_adj = cohort_AC/((cohort_max-cohort_na)/cohort_max)
            cohort_AC_per_case = cohort_AC_adj/len(cohort_cols)
            l_cohort_count.append(cohort_AC_per_case)
    df_svlen_bin["aou_count"] = l_aou_count
    df_svlen_bin["cohort_count"] = l_cohort_count
    
    #### Insertion
    df_svlen_bin['mid'] = (df_svlen_bin['start'] + df_svlen_bin['end']) / 2
    df_ins = df_svlen_bin[df_svlen_bin['mid'] > 0].copy()
    plt.figure(figsize=(6, 4), dpi=120)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(df_ins['mid'], df_ins['aou_count'], color='orange', label='Cohort', lw=1)
    plt.plot(df_ins['mid'], df_ins['cohort_count'], color='purple', lw=0.5, linestyle='--', label='All of Us')
    plt.grid(False)  # turn off all default grids
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    plt.xlim(sv_min,20000)
    plt.ylim(ymin=0.1)
    plt.xlabel('SV size (bp)')
    plt.ylabel('Count per case')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    output = os.path.join(outdir, 'svlen_INS.png')
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close() 
    
    #### Deletion
    df_del = df_svlen_bin[df_svlen_bin['mid'] < 0].copy()
    plt.figure(figsize=(6, 4), dpi=120)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(-df_del['mid'], df_del['aou_count'], color='orange', label='Cohort', lw=1)
    plt.plot(-df_del['mid'], df_del['cohort_count'], color='purple', lw=0.5, linestyle='--', label='All of Us')
    plt.xlim(20000, sv_min) 
    plt.ylim(ymin=0.1)
    plt.grid(False)
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=10))
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('SV size (bp)')
    plt.ylabel('Count per case')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    output = os.path.join(outdir, 'svlen_DEL.png')
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    
    ###################
    ## Impact of SVs on distinct genomic features
    ###################
    print("## Impact of SVs on distinct genomic features ##")
#     df_indels = df.loc[(df["SVTYPE"]=='INS') | (df["SVTYPE"]=='DEL')]
    df_indels = df.copy()
    df_indels["impact"] = ''

    dict_assign={'feature_elongation&coding_sequence_variant':'CDS',
    'feature_truncation&coding_sequence_variant&intron_variant':'CDS',
    'feature_truncation&inframe_deletion':'CDS',
    'stop_lost&feature_truncation&coding_sequence_variant&3_prime_UTR_variant&intron_variant':'CDS',
    'feature_truncation&coding_sequence_variant&5_prime_UTR_variant&intron_variant':'CDS',
    'splice_polypyrimidine_tract_variant&intron_variant':'CDS',
    'frameshift_variant&feature_truncation':'CDS',
    'feature_elongation&coding_sequence_variant&3_prime_UTR_variant':'CDS',
    'feature_truncation&coding_sequence_variant&5_prime_UTR_variant':'CDS',
    'splice_polypyrimidine_tract_variant&intron_variant&non_coding_transcript_variant':'CDS',
    'stop_lost&feature_truncation&coding_sequence_variant&5_prime_UTR_variant&3_prime_UTR_variant&intron_variant':'CDS',
    'stop_lost&feature_truncation&coding_sequence_variant&3_prime_UTR_variant':'CDS',
    'frameshift_variant&stop_lost&feature_truncation&3_prime_UTR_variant':'CDS',
    'stop_lost&feature_truncation&coding_sequence_variant&5_prime_UTR_variant&3_prime_UTR_variant':'CDS',
    'feature_elongation&coding_sequence_variant&5_prime_UTR_variant&3_prime_UTR_variant&intron_variant':'CDS',
    'frameshift_variant&start_lost&feature_truncation&start_retained_variant&5_prime_UTR_variant':'CDS',
    'start_lost&feature_elongation&start_retained_variant&coding_sequence_variant&intron_variant':'CDS',
    'stop_lost&feature_truncation&inframe_deletion':'CDS',
    'frameshift_variant&stop_lost&start_lost&feature_truncation&start_retained_variant&5_prime_UTR_variant&3_prime_UTR_variant':'CDS',
    'feature_truncation&coding_sequence_variant':'CDS',
    'feature_elongation&3_prime_UTR_variant':'UTR',
    'feature_truncation&3_prime_UTR_variant':'UTR',
    '3_prime_UTR_variant':'UTR',
    '5_prime_UTR_variant&intron_variant':'UTR',
    '5_prime_UTR_variant':'UTR',
    'feature_elongation&coding_sequence_variant&5_prime_UTR_variant&intron_variant':'UTR',
    'feature_elongation&coding_sequence_variant&3_prime_UTR_variant&intron_variant':'UTR',
    'coding_sequence_variant&3_prime_UTR_variant':'UTR',
    'feature_elongation&5_prime_UTR_variant':'UTR',
    'feature_truncation&5_prime_UTR_variant&intron_variant':'UTR',
    'coding_sequence_variant&5_prime_UTR_variant':'UTR',
    'feature_elongation&5_prime_UTR_variant&intron_variant':'UTR',
    'feature_truncation&5_prime_UTR_variant':'UTR',
    'coding_sequence_variant&5_prime_UTR_variant&3_prime_UTR_variant&intron_variant':'UTR',
    'start_lost&feature_elongation&start_retained_variant&coding_sequence_variant&5_prime_UTR_variant':'UTR',
    'feature_elongation&3_prime_UTR_variant&intron_variant':'UTR',
    'feature_truncation&3_prime_UTR_variant&intron_variant':'UTR',
    '3_prime_UTR_variant&intron_variant':'UTR',
    'intron_variant':'Intron',
    'feature_elongation&coding_sequence_variant&intron_variant':'Intron',
    'coding_sequence_variant&3_prime_UTR_variant&intron_variant':'Intron',
    'coding_sequence_variant&5_prime_UTR_variant&intron_variant':'Intron',
    'intron_variant&NMD_transcript_variant':'Intron',
    'regulatory_region_variant':'Intergenic', # TFBS
    'regulatory_region_ablation&regulatory_region_variant':'Intergenic', # TFBS
    'intron_variant&non_coding_transcript_variant':'non-coding gene',
    'feature_elongation&non_coding_transcript_exon_variant':'non-coding gene',
    'feature_truncation&non_coding_transcript_exon_variant&intron_variant':'non-coding gene',
    'feature_elongation&non_coding_transcript_exon_variant&intron_variant':'non-coding gene',
    'feature_truncation&non_coding_transcript_exon_variant':'non-coding gene',
    'non_coding_transcript_exon_variant&intron_variant':'non-coding gene',
    'non_coding_transcript_exon_variant':'non-coding gene',
    'non_coding_transcript_variant':'non-coding gene',
    'mature_miRNA_variant':'non-coding gene',
    'upstream_gene_variant':'Intergenic',
    'downstream_gene_variant':'Intergenic',
    'intergenic_variant':'Intergenic'}

    for key in dict_assign.keys():
        df_indels.loc[df_indels['Consequence']==key,'impact'] = dict_assign[key]

    _list=['protein_coding','IG_V_gene','IG_D_gene','TR_V_gene'
        ,'protein_coding_LoF','TR_C_gene','IG_C_gene']
    mask = (df_indels['Consequence']=='transcript_ablation')&(df_indels['BIOTYPE'].apply(lambda x: x in _list))
    df_indels.loc[mask,'impact'] = 'CDS'

    mask = (df_indels['Consequence']=='transcript_ablation')&(df_indels['BIOTYPE'].apply(lambda x: x not in _list))
    df_indels.loc[mask,'impact'] = 'non-coding gene'
    print(df_indels['impact'].value_counts())
    output = os.path.join(outdir,'')
    df_indels.to_csv(vep_csv.replace('.csv','.impact_annotated.csv'), index=False)
    
    ###################
    ## Plot SV Frequency: (1) Absolute count
    ###################
    print("## Plot SV Frequency: (1) Absolute count ##")    
    AC_aou = (df_indels[aou_cols]>=1).sum(axis=1)
    base_aou = (~df_indels[aou_cols].isna()).sum(axis=1)
    AC_cohort = (df_indels[cohort_cols]>=1).sum(axis=1)
    base_cohort = (~df_indels[cohort_cols].isna()).sum(axis=1)
    AC_all = (df_indels[aou_cols+cohort_cols]>=1).sum(axis=1)
    base_all = (~df_indels[aou_cols+cohort_cols].isna()).sum(axis=1)
    df_indels["af_aou"] = AC_aou/base_aou ## Allele frequecy
    df_indels["ac_aou"] = AC_aou ## Allele count
    df_indels["af_cohort"] = AC_cohort/base_cohort
    df_indels["ac_cohort"] = AC_cohort
    df_indels["af_all"] = AC_all/base_all
    df_indels["ac_all"] = AC_all
    VAF_cutoff = 0.01
    df_indels.loc[(df_indels['af_all']>=VAF_cutoff), 'category_all']='common'
    df_indels.loc[(df_indels['af_all']<VAF_cutoff), 'category_all']='rare'
    df_indels.loc[df_indels['ac_all']==1, 'category_all']='singleton'

    df_indels["_count"]=1
    pv_table = pd.pivot_table(df_indels[['impact','category_all','_count']], index='impact', columns='category_all', aggfunc = 'count')
    pv_table.columns = [c[1] for c in pv_table.columns]
    pv_table_cols = [col for col in ['Intergenic', 'Intron', 'non-coding gene', 'UTR', 'CDS'] if col in pv_table.columns.tolist()]
    pv_table = pv_table.loc[pv_table_cols]
    
    #### Absolute count
    fig, ax = plt.subplots(figsize=(8, 4), dpi=200)
    colors = ['#44722f', '#77a843', '#c5e0b4']
    bar_height = 0.7
    
    left_offset = np.zeros(len(pv_table))
    for i, col in enumerate(pv_table.columns): # Loop through columns to create stacked bars
        ax.barh(pv_table.index, pv_table[col], height=bar_height, left=left_offset, color=colors[i], label=col, edgecolor='white')
        left_offset += pv_table[col]
    ax.tick_params(axis='y', length=0, labelsize=14) # No y-axis ticks
    ax.tick_params(axis='x', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    total_counts = pv_table.sum(axis=1) # Add Total Count Labels 
    for i, total in enumerate(total_counts): # Loop through each bar to add the text
        ax.text(
            total + 50, # X-position: end of the bar + a small offset
            i, # Y-position: the index of the bar
            f' {total}', # The text to display (the total count)
            ha='left', # Horizontal alignment
            va='center', # Vertical alignment
            fontsize=12
        )
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.legend()
    output = os.path.join(outdir, 'SV_Frequency.png')
    plt.savefig(output, dpi=200, bbox_inches="tight")
    plt.close()
    
    ###################
    ## Plot SV Frequency: (2) Count in ratio
    ###################
    print("## Plot SV Frequency: (2) Count in ratio ##")    
    pv_table_ratio = pv_table.copy()
    pv_table_ratio['sum'] = pv_table_ratio.sum(axis=1)
    for col in pv_table_ratio.columns: pv_table_ratio[col] = pv_table_ratio[col]/pv_table_ratio['sum']
    fig, ax = plt.subplots(figsize=(8, 4), dpi=200)
    colors = ['#44722f', '#77a843', '#c5e0b4']
    bar_height = 0.7
    left_offset = np.zeros(len(pv_table_ratio)) # Loop through columns to create stacked bars
    for i, col in enumerate(pv_table.columns):
        ax.barh(pv_table_ratio.index, pv_table_ratio[col], height=bar_height, left=left_offset, color=colors[i], label=col, edgecolor='white')
        left_offset += pv_table_ratio[col]
    ax.tick_params(axis='y', length=0, labelsize=14) # No y-axis ticks
    ax.tick_params(axis='x', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    output = os.path.join(outdir, 'SV_Frequency(ratio).png')
    plt.savefig(output, dpi=200, bbox_inches="tight")
    plt.close()
    
    ###################
    ## Aou vs Cohort: By Impact x Frequency category
    ###################
    print("## Aou vs Cohort: By Impact x Frequency category ##")    
    df_indels_aou = df_indels.loc[(df_indels["ac_aou"]>=1),['impact','category_all']]
    df_indels_cohort = df_indels.loc[(df_indels["ac_cohort"]>=1),['impact','category_all']]
    df_indels_aou["_count"]=1
    df_indels_cohort["_count"]=1

    pv_table_aou = pd.pivot_table(df_indels_aou, index='impact', columns='category_all', aggfunc = 'count')
    pv_table_aou.columns = [c[1] for c in pv_table_aou.columns]
    pv_table_aou = pv_table_aou.loc[pv_table_cols]

    pv_table_cohort = pd.pivot_table(df_indels_cohort, index='impact', columns='category_all', aggfunc = 'count')
    pv_table_cohort.columns = [c[1] for c in pv_table_cohort.columns]
    pv_table_cohort = pv_table_cohort.loc[pv_table_cols]

    l_impact = []
    l_cat =[]
    l_ratio = []
    l_base = []
    l_count = []
    l_group = []

    for impact in pv_table_cols:
        for cat in ['common','rare','singleton']:
            for tab, group in zip([pv_table_aou, pv_table_cohort],['AoU','FPC']):
                l_impact.append(impact)
                l_cat.append(cat)
                l_group.append(group)

                ## calc ratio
                base = tab[cat].sum()
                count = tab.loc[impact,cat].sum()
                ratio = count/base

                l_ratio.append(ratio)
                l_base.append(base)
                l_count.append(count)
    df_summary=pd.DataFrame({
         "impact":l_impact
        ,"cat":l_cat
        ,"ratio":l_ratio
        ,"base":l_base
        ,"count":l_count
        ,"group":l_group
    }) 


    for impact in pv_table_cols:
        plt.figure(figsize=(3,3),dpi=200)

        df_draw = df_summary.loc[df_summary["impact"]==impact,['cat','ratio','group']]
        df_draw = pd.pivot(df_draw,columns="group",index="cat")
        df_draw = df_draw.reset_index()
        df_draw.columns = ['category']+[c[1] for c in df_draw.columns[1:]]

        fig, ax = plt.subplots(figsize=(4, 4), dpi=200) # Increased DPI for publication quality
        bar_width = 0.35
        x_pos = np.arange(len(df_draw['category']))
        color_series1 = 'orange'
        color_series2 = 'purple'

        ax.bar(x_pos - bar_width*1.1/2, df_draw['AoU']*100, bar_width, label='AoU', color=color_series1)
        ax.bar(x_pos + bar_width*1.1/2, df_draw['FPC']*100, bar_width, label='FPC', color=color_series2)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(df_draw['category'], fontsize=14) # Larger font for readability
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_edgecolor('k') # Light gray for the left spine
        ax.spines['bottom'].set_edgecolor('k') # Light gray for the bottom spine
        ax.grid(axis='y', linestyle='-', alpha=0.6, color='lightgray')
        ax.set_axisbelow(True) # Ensures grid lines are behind the bars
        ax.tick_params(axis='y', labelsize=12, colors='k') # Darker, smaller ticks
        ax.tick_params(axis='x', length=0) # Remove x-axis tick marks
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=12)
        plt.tight_layout()
        plt.title(impact)
        output = os.path.join(outdir, 'SV_ratio_by_region_by_freq_by_group__{}.png'.format(impact))
        plt.savefig(output, dpi=200, bbox_inches="tight")
        plt.close()
    
    ###################
    ## Chi-sq test
    ###################
    print("## Chi-sq test ##")    
    print()
    print("impact\tcategory\tp_value")
    for impact in pv_table_cols:
        for cat in ["common","rare","singleton"]:
            mask1 = df_summary["impact"]==impact
            mask2 = df_summary["cat"]==cat

            x_aou = df_summary.loc[mask1&mask2&(df_summary["group"]=="AoU")]
            xx_aou = [x_aou["count"].iloc[0], x_aou["base"].iloc[0] - x_aou["count"].iloc[0]]

            x_fpc = df_summary.loc[mask1&mask2&(df_summary["group"]=="FPC")]
            xx_fpc = [x_fpc["count"].iloc[0], x_fpc["base"].iloc[0] - x_fpc["count"].iloc[0]]

            observed_data = np.array([xx_aou, xx_fpc])
            chi2_stat, p_value, dof, expected_values = chi2_contingency(observed_data)
            p_value = round(p_value,3)
            print(f"{impact}\t{cat}\t{p_value}")


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Post VEP analysis 1")
    parser.add_argument('--vep_csv', help="parsed vep in csv")
    parser.add_argument('--precede_mastersheet', help="precede master spreadsheet")
    parser.add_argument('--outdir', help="output directory")
    args = parser.parse_args()
    
    main(args.vep_csv, args.precede_mastersheet, args.outdir)

    

    

