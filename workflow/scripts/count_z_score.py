import numpy as np
import pandas as pd


def contig_stats(group):
    total_weight = group['length'].sum()
    weighted_mean_cov = np.average(group['cov'], weights=group['length'])
    weighted_std_cov = np.sqrt(np.average((group['cov'] - weighted_mean_cov)**2, weights=group['length']))
    return pd.Series({
        'contig_length': total_weight,
        'weighted_mean_cov': weighted_mean_cov,
        'weighted_std_cov': weighted_std_cov
    })
    
def merge_z_score_counting(interval_coverage_df, contig_stats_df):
    # Merge the two DataFrames on the 'contig_header' column
    merged_df = pd.merge(interval_coverage_df, contig_stats_df, on='contig_header')
    # Calculate the z-score for each interval
    merged_df['z_score'] = (merged_df['cov'] - merged_df['weighted_mean_cov']) / merged_df['weighted_std_cov']
    # Handle division by zero if needed
    merged_df.loc[merged_df['weighted_std_cov'] == 0, 'z_score'] = 0
    # Drop unnecessary columns
    columns_to_drop = ['weighted_mean_cov', 'weighted_std_cov', 'contig_length']
    merged_df.drop(columns=columns_to_drop, inplace=True)
    return merged_df