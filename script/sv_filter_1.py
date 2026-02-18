import pandas as pd
import numpy as np
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool
import os
import gzip
import argparse

def process_gt_column(col_data_tuple):
    """
    Worker function to process a single column.
    
    Args:
        col_data_tuple (tuple): A tuple containing the column name (str) and
                                the column data (pd.Series).

    Returns:
        tuple: A tuple containing the column name and the processed Series.
    """
    col_name, data, idx = col_data_tuple
    gt_data = data.str.split(':').str[idx]
    return col_name, gt_data

def create_df_gt_multiprocess(df, sample_cols, num_threads, idx):
    """
    Creates the df_gt DataFrame using multiprocessing.
    
    Args:
        df (pd.DataFrame): The input DataFrame.
        sample_cols (list): A list of column names to process.
        num_threads (int): The number of threads to use.
    
    Returns:
        pd.DataFrame: The resulting df_gt DataFrame.
    """
    # Prepare input for the pool by creating a list of (column_name, Series) tuples
    column_data = [(col, df[col], idx) for col in sample_cols]

    # Use a multiprocessing pool to parallelize the processing
    with Pool(num_threads) as p:
        # The tqdm wrapper provides a progress bar
        results = list(tqdm(p.imap(process_gt_column, column_data), total=len(column_data), desc="Processing columns"))

    # Reconstruct the dictionary from the parallel results
    _df_gt = {col_name: series for col_name, series in results}

    # Create the final DataFrame and ensure columns are in the correct order
    df_gt = pd.DataFrame(_df_gt)[sample_cols]
    
    return df_gt

def variant_filter(vcf, case, ctrl, outdir, num_threads):
    filt_variant_bed = os.path.join(outdir,'sv_filter1.bed')
    filt_variant_id = os.path.join(outdir,'sv_filter1.id')
    filt_sample =  os.path.join(outdir,'filt_sample.txt')
    
    
    ###################
    # Read VCF
    ###################
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
    df_vcf = df_vcf.rename(columns={'#CHROM': 'CHROM'})

    df_ctrl =pd.read_csv(ctrl)
    df_case =pd.read_csv(case)
    sample_case = [sample for sample in  df_case['analysis_sample_id'].astype(str) if sample in df_vcf.columns.tolist()]
    sample_ctrl = [f'{sample}.sniffles2.hg38' for sample in  df_ctrl['analysis_sample_id'].astype(str) if f'{sample}.sniffles2.hg38' in df_vcf.columns.tolist()]
    sample_all= sample_case + sample_ctrl 

    def parse_info(info, keyword):
        if keyword in info: return info.split(keyword+'=')[1].split(';')[0]
        else: return np.nan

    for kw,dtype in zip(['AC','SVLEN','SVTYPE','STRAND','SUPPORT'],['float','float','str','str','float']):
            df_vcf[kw] = df_vcf["INFO"].apply(lambda x: parse_info(x, kw)).astype(dtype)
    
    ###################
    # Filtering
    ###################
    print('Before filtering (Num variants):')
    print(df_vcf.shape[0])

    print('Mask1: GT missing rate < 20% at locus After setting DP<8 or GQ<10 as missing value')
    df_gt = create_df_gt_multiprocess(df_vcf, sample_all, num_threads, 0)
    df_gq = create_df_gt_multiprocess(df_vcf, sample_all, num_threads, 1)
    df_dr = create_df_gt_multiprocess(df_vcf, sample_all, num_threads, 2)
    df_dv = create_df_gt_multiprocess(df_vcf, sample_all, num_threads, 3)

    d = {'./.':np.nan, '0/0':0, '0/1':1, '1/1':2}
    for col in sample_all:
        df_gt[col] = df_gt[col].map(d)

    df_gq = df_gq.astype(int)
    df_dr = df_dr.astype(int)
    df_dv = df_dv.astype(int)
    df_dp = df_dr+df_dv

    df_gt[df_dp<8]=np.nan
    mask_gq_1 = df_gq<10
    mask_gq_2 = (df_gt==0)&(df_gq==0)
    mask_gq = mask_gq_1 & (~mask_gq_2)
    df_gt[mask_gq]=np.nan

    mask1 = (df_gt.isna().sum(axis=1)/len(sample_all))<0.2
    print(sum(mask1))

    print('Mask2: AC>=1 @ cohort_cols')
    df_vcf_case = df_vcf[sample_case].copy()
    mask2 = []
    for i in range(df_vcf_case.shape[0]):
        uniq = df_vcf_case.loc[i].str.split(':').str[0].unique()
        uniq = set(uniq)
        if len(uniq - set(['./.', '0/0']))>0: mask2.append(True)
        else: mask2.append(False)
    mask2 = np.array(mask2)
    print(sum(mask1&mask2))

    print('Mask3: SVLEN -- (df["SVLEN"].isna() | (df["SVLEN"].abs()>=50))')
    mask3 = (df_vcf["SVLEN"].isna() | (df_vcf["SVLEN"].abs()>=50))
    print(sum(mask1&mask2&mask3))

    print('Mask4: Unplaced Chr')
    mask4 = df_vcf["CHROM"].apply(lambda x: 'chrUn' not in x)& df_vcf["CHROM"].apply(lambda x: 'random' not in x)
    print(sum(mask1&mask2&mask3&mask4))
    
    ###################
    # Saving
    ###################
    df_vcf_filt = df_vcf.loc[mask1&mask2&mask3&mask4].reset_index(drop=True)
    df_vcf_filt_bed = df_vcf_filt[['CHROM','POS','POS']].copy()
    df_vcf_filt_bed.columns = ['CHROM','START','END']
    df_vcf_filt_bed['START'] = df_vcf_filt_bed['START']-1
    df_vcf_filt_bed.to_csv(filt_variant_bed, sep="\t", index=False, header=None)
    print('Saving: ', filt_variant_bed)
    print('saving: ', filt_sample)

    with open(filt_sample, 'w') as f:
        for s in sample_case: f.write(s+'\n')
        for s in sample_ctrl: f.write(s+'\n')
    with open(filt_variant_id, 'w') as f:
        for ID in df_vcf_filt["ID"]: f.write(ID+'\n')
            
            


#     print(df.shape[0])
#     mask1 = (df["AC"]>=1)&(df["SVLEN"].isna() | (df["SVLEN"].abs()>=50))
#     print(df.loc[mask1].shape[0])
# #     mask2 = df["STRAND"] == "+-"
# #     print(df.loc[mask1&mask2].shape[0])
#     mask3 = df["CHROM"].apply(lambda x: 'chrUn' not in x)& df["CHROM"].apply(lambda x: 'random' not in x)
#     print(df.loc[mask1&mask3].shape[0]) ## (-) mask2&

#     print(df.value_counts("SVTYPE"))
#     df_filt1 = df.loc[mask1&mask3].reset_index(drop=True) ## (-) mask2&
#     print(df_filt1.value_counts("SVTYPE"))
#     # Execute the multiprocessing function
#     df_gt = create_df_gt_multiprocess(df_filt1, sample_cols, num_threads, 0)
#     df_gq = create_df_gt_multiprocess(df_filt1, sample_cols, num_threads, 1)
#     df_dr = create_df_gt_multiprocess(df_filt1, sample_cols, num_threads, 2)
#     df_dv = create_df_gt_multiprocess(df_filt1, sample_cols, num_threads, 3)

#     d = {'./.':np.nan, '0/0':0, '0/1':1, '1/1':2}
#     for col in tqdm(sample_cols):
#         df_gt[col] = df_gt[col].map(d)

#     mask = (df_gt.isna().sum(axis=1)/len(sample_cols))<0.2
#     df_filt2 = df_filt1.loc[mask]
#     print(df_filt2.value_counts("SVTYPE"))
#     df_gt2 = df_gt.loc[mask]
#     df_gq2 = df_gq.loc[mask].astype(int)
#     df_dr2 = df_dr.loc[mask].astype(int)
#     df_dv2 = df_dv.loc[mask].astype(int)
#     df_dp2 = df_dr2+df_dv2

#     print("Num variant Loci: Before ({}) - After({}) ; -{}%"\
#       .format(df.shape[0], df_gt2.shape[0], round(100*(1-df_gt2.shape[0]/df.shape[0]),2)))


#     ###################
#     # Filtering : Num usable samples (calc on rare variant loci)
#     ###################
#     # Define a cutoff for "rare" loci
#     mask1 = (df_gt2[aou_cols + cohort_cols] > 0).sum(axis=1) / len(aou_cols + cohort_cols) < 0.01
#     mask2 = (df_gt2[aou_cols + cohort_cols] > 0).sum(axis=1) / len(aou_cols + cohort_cols) > 0.99
#     rare_loci_mask = mask1 | mask2
#     num_rare_loci = sum(rare_loci_mask)

#     # Initialize lists to store ratios for different depth cutoffs
#     l_cutoff = []
#     l_aou_count = []
#     l_cohort_count = []

#     # Iterate through depth cutoffs
#     for cutoff in range(23):
#         l_cutoff.append(cutoff)
#         aou_covered_loci_counts = (df_dp2.loc[rare_loci_mask, aou_cols] >= cutoff).sum(axis=1)
#         l_aou_count.append(aou_covered_loci_counts.tolist())

#         cohort_covered_loci_counts = (df_dp2.loc[rare_loci_mask, cohort_cols] >= cutoff).sum(axis=1)
#         l_cohort_count.append(cohort_covered_loci_counts.tolist())

#     # ---
#     # Plotting
#     # ---

#     # Set a consistent style for the plot
#     sns.set_style("whitegrid")
#     plt.style.use('seaborn-v0_8-deep')

#     # Transform the data for plotting
#     df_aou = pd.DataFrame(l_aou_count, index=l_cutoff).T.melt(var_name='Depth Cutoff', value_name='Count')
#     df_aou['Cohort'] = 'All of Us'

#     df_cohort = pd.DataFrame(l_cohort_count, index=l_cutoff).T.melt(var_name='Depth Cutoff', value_name='Count')
#     df_cohort['Cohort'] = 'Cohort'

#     # Concatenate the dataframes for a unified plot
#     df_plot = pd.concat([df_aou, df_cohort], ignore_index=True)

#     # Create the plot
#     plt.figure(figsize=(4, 4), dpi=150)
#     ax = sns.lineplot(
#         data=df_plot.loc[df_plot['Depth Cutoff']>0],
#         x='Depth Cutoff',
#         y='Count',
#         hue='Cohort',
#         style='Cohort',
#         markers=True,
#         dashes=False,
#         markersize=8,
#         linewidth=2.5,
#     )

#     # Customize the plot for publication
#     ax.set_xlabel('Minimum Sequencing Depth', fontsize=12, labelpad=10)
#     ax.set_ylabel('Num of usuable samples', fontsize=12, labelpad=10)
#     ax.tick_params(labelsize=10)

#     # Set x-axis ticks to be cleaner
#     ax.set_xticks(range(0, 24, 2))

#     # Add a legend and position it appropriately
#     ax.legend(loc='lower left', frameon=True, fontsize=10, title_fontsize='12')

#     # Final touches for a clean look
#     sns.despine(trim=True)
#     plt.tight_layout()

#     # Save the figure to a file
#     plt.savefig(os.path.join(outdir,'depth_cutoff.png'), dpi=300, bbox_inches='tight')

#     # APPLY FILTER -----> DP>=8
#     df_gt3 = df_gt2.copy()
#     df_gt3[df_dp2<8]=np.nan
#     mask = (df_gt3.isna().sum(axis=1)/len(sample_cols))<0.2


#     df_filt3 = df_filt2.loc[mask].copy()
#     print(df_filt3.value_counts("SVTYPE"))
#     df_gt3 = df_gt3.loc[mask].copy()
#     df_gq3 = df_gq2.loc[mask].copy()
#     df_dr3 = df_dr2.loc[mask].copy()
#     df_dv3 = df_dv2.loc[mask].copy()
#     df_dp3 = df_dp2.loc[mask].copy()

#     print("Num variant loci: Before ({}) - After({}) ; -{}%"\
#           .format(df_gt2.shape[0], df_gt3.shape[0], round(100*(1-df_gt3.shape[0]/df_gt2.shape[0]),2)))


#     ###################
#     # Filtering: % ALT calls to be utilized at rare variant loci depending on GQ filter
#     ###################
#     # Define a cutoff for "rare" loci
#     mask1 = (df_gt3[aou_cols + cohort_cols] > 0).sum(axis=1) / len(aou_cols + cohort_cols) < 0.01
#     mask2 = (df_gt3[aou_cols + cohort_cols] > 0).sum(axis=1) / len(aou_cols + cohort_cols) > 0.99
#     rare_loci_mask = mask1 | mask2
#     num_rare_loci = sum(rare_loci_mask)

#     # Initialize lists to store ratios for different depth cutoffs
#     l_cutoff = []
#     l_aou_count = []
#     l_cohort_count = []

#     ALT = df_gt3.loc[rare_loci_mask, aou_cols]>=1
#     gq_values = df_gq3.loc[rare_loci_mask, aou_cols][ALT].values.ravel()
#     gq_values_aou = gq_values[np.isfinite(gq_values)]

#     ALT = df_gt3.loc[rare_loci_mask, cohort_cols]>=1
#     gq_values = df_gq3.loc[rare_loci_mask, cohort_cols][ALT].values.ravel()
#     gq_values_cohort = gq_values[np.isfinite(gq_values)]


#     # Iterate through depth cutoffs (e.g., 0 to 19)
#     for cutoff in tqdm(range(60)):
#         l_cutoff.append(cutoff)

#         n_alt_calls_aou = sum(gq_values_aou>=cutoff)/len(gq_values_aou)
#         l_aou_count.append(n_alt_calls_aou)

#         n_alt_calls_cohort = sum(gq_values_cohort>=cutoff)/len(gq_values_cohort)
#         l_cohort_count.append(n_alt_calls_cohort)
#     # ---
#     # Plotting
#     # ---

#     # Set a consistent style for the plot
#     sns.set_style("whitegrid")
#     plt.style.use('seaborn-v0_8-deep')

#     # Transform the data for plotting
#     df_aou = pd.DataFrame(l_aou_count, index=l_cutoff).T.melt(var_name='GQ Cutoff', value_name='Count')
#     df_aou['Cohort'] = 'All of Us'

#     df_cohort = pd.DataFrame(l_cohort_count, index=l_cutoff).T.melt(var_name='GQ Cutoff', value_name='Count')
#     df_cohort['Cohort'] = 'Cohort'

#     # Concatenate the dataframes for a unified plot
#     df_plot = pd.concat([df_aou, df_cohort], ignore_index=True)

#     # Create the plot
#     plt.figure(figsize=(10, 4), dpi=150)
#     ax = sns.lineplot(
#         data=df_plot,
#         x='GQ Cutoff',
#         y='Count',
#         hue='Cohort',
#         style='Cohort',
#         markers=True,
#         dashes=False,
#         markersize=8,
#         linewidth=2.5,
#     )

#     # Customize the plot for publication
#     ax.set_xlabel('Minimum GQ', fontsize=12, labelpad=10)
#     ax.set_ylabel('Num of ALT calls', fontsize=12, labelpad=10)
#     ax.tick_params(labelsize=10)

#     # Set x-axis ticks to be cleaner
#     ax.set_xticks(range(0, 60, 2))

#     ax.set_ylim(0,1.03)

#     # Final touches for a clean look
#     sns.despine(trim=True)
#     plt.tight_layout()

#     #Save the figure to a file
#     plt.savefig(os.path.join(outdir,'GQ_cutoff.png'), dpi=300, bbox_inches='tight')

#     df_gt4 = df_gt3.copy()
#     mask1 = df_gq3<20
#     mask2 = (df_gt3==0)&(df_gq3==0)
#     mask = mask1 & (~mask2)
#     df_gt4[mask]=np.nan

#     mask = (df_gt4.isna().sum(axis=1)/len(sample_cols))<0.2
#     mask = mask & (df_gt4.sum(axis=1)>0)

#     df_filt4 = df_filt3.loc[mask].copy()
#     print(df_filt4.value_counts("SVTYPE"))
#     df_gt4 = df_gt4.loc[mask].copy()
#     df_gq4 = df_gq3.loc[mask].copy()
#     df_dr4 = df_dr3.loc[mask].copy()
#     df_dv4 = df_dv3.loc[mask].copy()
#     df_dp4 = df_dp3.loc[mask].copy()

#     print("Num variant loci: Before ({}) - After({}) ; -{}%"\
#           .format(df_gt3.shape[0], df_gt4.shape[0], round(100*(1-df_gt4.shape[0]/df_gt3.shape[0]),2)))

#     pd.concat((df_filt4[["CHROM","POS","ID"]], ~df_gt4.isna()),axis=1).to_csv(outcsv, index=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Filter VCF")
    parser.add_argument('--vcf', help="vep result in vcf format")
    parser.add_argument('--case')
    parser.add_argument('--ctrl')
    parser.add_argument('--outdir', help="outdir for figures")
    parser.add_argument('--num_threads')
    args = parser.parse_args()
    
    variant_filter(args.vcf, args.case, args.ctrl, args.outdir, int(args.num_threads))
