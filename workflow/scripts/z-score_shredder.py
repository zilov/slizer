#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#@created: 13.10.2023
#@author: Danil Zilov
#@contact: zilov.d@gmail.com

import pandas as pd
import argparse
from count_z_score import contig_stats, merge_z_score_counting
import logging

def setup_logging(output_folder, prefix):
    logging.basicConfig(filename=f"{output_folder}/{prefix}_analysis.log",
                        level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] - %(message)s")


def fasta_reader_yield(path_to_fasta_file):
    """
    Reads a BED file containing contig coverage information and calculates z-scores.
    
    Parameters:
        path_to_fasta_file (str): Path to the FASTA file.
        
    Returns:
        generator which returns header: sequence
    """
    header = None
    with open(path_to_fasta_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header,"".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header,"".join(seq)


def read_bed_file(bed_file_path):
    logging.info("Read contigs coverage...")
    icov = pd.read_csv(bed_file_path, sep='\t', header=None, names=['contig_header', 'start', 'end', 'cov'])
    icov['length'] = icov['end'] - icov['start']
    return icov


# Function to calculate statistics for each contig
def calculate_contig_stats(icov):
    logging.info("Counting contig stats...")
    contig_stats_df = icov.groupby('contig_header').apply(contig_stats).reset_index()
    return contig_stats_df


# Function to merge z-score to interval DataFrame
def merge_z_score_to_intervals(icov, contig_stats_df):
    logging.info("Merging z-scores...")
    interval_cov_with_zscore_df = merge_z_score_counting(icov.copy(), contig_stats_df.copy())
    return interval_cov_with_zscore_df


# Function to calculate the global mean coverage
def calculate_global_mean_coverage(icov):
    global_mean_coverage = (icov['cov'] * icov['length']).sum() / icov['length'].sum()
    return global_mean_coverage


def filter_low_coverage_regions(interval_df, length_threshold, tail_threshold, z_threshold, cov_threshold):
    """
    Filters out low coverage regions based on given thresholds for length, tail position, and z-score.
    Adds an extra condition to not remove regions if their coverage is greater than a global mean coverage threshold.

    Parameters:
        interval_df (pd.DataFrame): DataFrame containing intervals and their z-scores.
        length_threshold (int): Minimum length for low-coverage regions.
        tail_threshold (int): Size of the tail regions to exclude.
        z_threshold (float): Z-score threshold for low coverage.
        cov_threshold (float): Coverage threshold based on global mean.

    Returns:
        pd.DataFrame: DataFrame containing only the low-coverage regions that meet the criteria.
    """
    
    # Filter by length
    length_filtered_df = interval_df[interval_df['length'] > length_threshold]
    
    # Filter by global mean coverage
    mean_cov_filtered_df = length_filtered_df[length_filtered_df['cov'] < cov_threshold]

    # Filter by z-score
    z_filtered_df = mean_cov_filtered_df[mean_cov_filtered_df['z_score'] < z_threshold]

    # Remove regions in the tails of the contigs
    final_filtered_df = z_filtered_df.copy()
    for contig in z_filtered_df['contig_header'].unique():
        contig_intervals = z_filtered_df[z_filtered_df['contig_header'] == contig]
        contig_length = contig_intervals['end'].max()  # Assuming the intervals cover the whole contig

        # Remove tail regions
        final_filtered_df = final_filtered_df[
            ~((final_filtered_df['contig_header'] == contig) & 
              ((final_filtered_df['start'] < tail_threshold) | 
               (final_filtered_df['end'] > contig_length - tail_threshold)))
        ]

    return final_filtered_df


def indicate_high_coverage_regions(interval_df, z_threshold):
    """
    Identifies high coverage regions based on a given z-score threshold.

    Parameters:
        interval_df (pd.DataFrame): DataFrame containing intervals and their z-scores.
        z_threshold (float): Z-score threshold for high coverage.

    Returns:
        pd.DataFrame: DataFrame containing only the high-coverage regions that meet the z-score criteria.
    """
    
    # Filter by z-score
    high_cov_regions = interval_df[interval_df['z_score'] > z_threshold]
    
    return high_cov_regions


def extract_low_coverage_sequences(fasta_reader, low_cov_regions):
    """
    Extract sequences corresponding to low-coverage regions and save to a dictionary.

    Parameters:
        fasta_reader (generator): Generator that yields tuples of (header, sequence) from a FASTA file.
        low_cov_regions (pd.DataFrame): DataFrame containing only the low-coverage regions.

    Returns:
        dict: Dictionary containing the extracted low-coverage sequences.
    """

    low_cov_fasta_dict = {}
    for record_id, sequence in fasta_reader:
        low_coverage_intervals = low_cov_regions[low_cov_regions['contig_header'] == record_id]
        
        # If there are no low-coverage regions for this contig, continue to the next iteration
        if low_coverage_intervals.empty:
            continue
        
        split_index = 1
        for _, row in low_coverage_intervals.iterrows():
            # Extract the sequence for the low-coverage region
            low_cov_seq = sequence[row['start']:row['end']]
            
            # Save to the dictionary
            low_cov_fasta_dict[f"{record_id}_lowcov_{split_index}_{row['start']}_{row['end']}_{row['cov']}_{row['z_score']}"] = low_cov_seq
            
            split_index += 1

    return low_cov_fasta_dict


def extract_splitted_fasta(fasta_reader, low_cov_regions, remove_low_coverage=False):
    """
    Extract sequences with or without low-coverage regions based on the 'hardcut' parameter and save to a dictionary.

    Parameters:
        fasta_reader (generator): Generator that yields tuples of (header, sequence) from a FASTA file.
        low_cov_regions (pd.DataFrame): DataFrame containing only the low-coverage regions.
        hardcut (bool): If True, completely remove low-coverage regions. If False, split at the end value of low-coverage regions.

    Returns:
        dict: Dictionary containing the extracted sequences with modifications based on 'hardcut' parameter.
    """
    
    splitted_fasta_dict = {}
    for record_id, sequence in fasta_reader:
        low_coverage_intervals = low_cov_regions[low_cov_regions['contig_header'] == record_id]
        
        if low_coverage_intervals.empty:
            splitted_fasta_dict[record_id] = sequence
            continue
        
        last_end = 0
        split_index = 1
        for _, row in low_coverage_intervals.iterrows():
            if remove_low_coverage:
                new_seq = sequence[last_end:row['start']]
            else:
                # Include the low-coverage region in the sequence if not hardcutting
                new_seq = sequence[last_end:row['end']]

            new_record_id = f"{record_id}_{split_index}"
            splitted_fasta_dict[new_record_id] = new_seq
            last_end = row['end']
            split_index += 1

        # Save the remaining part of the contig after the last low-coverage region
        if last_end < len(sequence):  # Ensure there's a remaining part to be saved
            new_seq = sequence[last_end:]
            new_record_id = f"{record_id}_{split_index}"
            splitted_fasta_dict[new_record_id] = new_seq

    return splitted_fasta_dict


def write_fasta_file(fasta_dict, out_file):
    with open(out_file, 'w') as fw:
        for h, s in fasta_dict.items():
            fw.write(f">{h}\n{s}\n")


def write_tables(low_cov_regions, high_cov_regions, contig_stats, low_cov_file, high_cov_file, contig_stats_file):
    """
    Write tables with low and high coverage region details and contig statistics to files.

    Parameters:
        low_cov_regions (pd.DataFrame): DataFrame containing only the low-coverage regions.
        high_cov_regions (pd.DataFrame): DataFrame containing only the high-coverage regions.
        contig_stats (pd.DataFrame): DataFrame containing statistics for each contig.
        low_cov_file (str): Path to the file where low-coverage region details will be saved.
        high_cov_file (str): Path to the file where high-coverage region details will be saved.
        contig_stats_file (str): Path to the file where contig statistics will be saved.
    """
    
    # Write low-coverage regions to file
    low_cov_regions.to_csv(low_cov_file, sep='\t', index=False)
    
    # Write high-coverage regions to file
    high_cov_regions.to_csv(high_cov_file, sep='\t', index=False)
    
    # Write contig statistics to file
    contig_stats.to_csv(contig_stats_file, sep='\t', index=False)


def write_final_statistics_report(low_cov_regions, high_cov_regions, output_file):
    """
    Write final statistics report including the total number of low/high-coverage regions,
    total removed sequences length, and total contigs splitted.

    Parameters:
        low_cov_regions (pd.DataFrame): DataFrame containing only the low-coverage regions.
        high_cov_regions (pd.DataFrame): DataFrame containing only the high-coverage regions.
        output_file (str): Path to the file where the final statistics will be saved.
    """
    
    total_low_cov_regions = len(low_cov_regions)
    total_high_cov_regions = len(high_cov_regions)
    total_removed_sequences_length = low_cov_regions['length'].sum()
    total_contigs_splitted = len(low_cov_regions['contig_header'].unique())
    
    with open(output_file, 'w') as f:
        logging.info(f"Total low-coverage regions:{total_low_cov_regions}")
        f.write(f"Total low-coverage regions:\t{total_low_cov_regions}\n")
        logging.info(f"Total high-coverage regions: {total_high_cov_regions}")
        f.write(f"Total high-coverage regions:\t{total_high_cov_regions}\n")
        logging.info(f"Total removed sequences length: {total_removed_sequences_length}")
        f.write(f"Total removed sequences length:\t{total_removed_sequences_length}\n")
        logging.info(f"Total contigs splitted: {total_contigs_splitted}") 
        f.write(f"Total contigs splitted:\t{total_contigs_splitted}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze and process genome assembly based on coverage data.")
    parser.add_argument("-b", "--bed_file", required=True, help="Path to the BED file containing interval coverage data.")
    parser.add_argument("-a", "--fasta_file", required=True, help="Path to the FASTA file containing the genome assembly.")
    parser.add_argument("-o", "--output_folder", required=True, help="Folder where all output files will be saved.")
    parser.add_argument("-p", "--prefix", default="", help="Prefix for all output files.")
    parser.add_argument("-c", "--mean_cov_fraction", type=float, default=0.25, 
                        help="Fraction of the global mean coverage to use as a threshold for filtering intervals."
                             "Intervals with coverage below fraction will be used for further analysis.")
    parser.add_argument("-l", "--length_threshold", type=int, default=500, help="Minimum length for low-coverage regions.")
    parser.add_argument("-t", "--tail_threshold", type=int, default=1000, help="Size of the tail regions to exclude.")
    parser.add_argument("-zl", "--z_threshold_low", type=float, default=-2, help="Z-score threshold for low coverage.")
    parser.add_argument("-zh", "--z_threshold_high", type=float, default=2, help="Z-score threshold for high coverage.")
    parser.add_argument("-r", "--remove_low_coverage", type=bool, default=False, help="Remove low coverage regions from splitted assembly/")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    # Fetch arguments
    bed_file_path = args.bed_file
    fasta_file_path = args.fasta_file
    output_folder = args.output_folder
    prefix = args.prefix
    length_threshold = args.length_threshold
    tail_threshold = args.tail_threshold
    z_threshold_low = args.z_threshold_low
    z_threshold_high = args.z_threshold_high
    mean_cov_fraction = args.mean_cov_fraction
    remove_low_coverage = args.remove_low_coverage
    
    # Output files
    low_cov_file = f"{output_folder}/{prefix}_low_cov_regions.tsv"
    high_cov_file = f"{output_folder}/{prefix}_high_cov_regions.tsv"
    contig_stats_file = f"{output_folder}/{prefix}_contig_stats.tsv"
    final_stats_file = f"{output_folder}/{prefix}_final_statistics_report.tsv"
    low_cov_fasta_out = f"{output_folder}/{prefix}_low_cov_regions.fa"
    splitted_fasta_out = f"{output_folder}/{prefix}_splitted.fa"
    
    # Initialize logging
    setup_logging(output_folder, prefix)
    
    logging.info("Starting analysis.")
    logging.info(f"""Run parameters are:
        bed_file_path: {bed_file_path},
        fasta_file_path = {fasta_file_path}
        output_folder = {output_folder}
        prefix = {prefix}
        length_threshold = {length_threshold}
        tail_threshold = {tail_threshold}
        z_threshold_low = {z_threshold_low}
        z_threshold_high = {z_threshold_high}
        mean_cov_fraction = {mean_cov_fraction}""")
    
    
    # Step 1: Read BED file and calculate Z-score
    logging.info("Reading BED file and calculating Z-score and mean coverage.")
    icov = read_bed_file(bed_file_path)
    contig_stats_df = calculate_contig_stats(icov)
    interval_cov_with_zscore_df = merge_z_score_to_intervals(icov, contig_stats_df)
    global_mean_coverage = calculate_global_mean_coverage(icov)
    logging.info(f"Mean coverage is {global_mean_coverage}")
    coverage_threshold = global_mean_coverage * mean_cov_fraction
    logging.info(f"Coverage threshold is {coverage_threshold}")

    # Step 2: Filter low coverage regions
    logging.info("Filtering low coverage regions.")
    low_cov_regions = filter_low_coverage_regions(interval_cov_with_zscore_df, length_threshold, tail_threshold, z_threshold_low, coverage_threshold)
    logging.info(f"Low coverage intervals found: {len(low_cov_regions)}")

    # Step 3: Indicate high coverage regions
    logging.info("Indicating high coverage regions.")
    high_cov_regions = indicate_high_coverage_regions(interval_cov_with_zscore_df, z_threshold_high)
    logging.info(f"High coverage intervals found: {len(low_cov_regions)}")
    
    # Step 4: Extract low coverage sequences and save to a dictionary
    fasta_reader = fasta_reader_yield(fasta_file_path)
    low_cov_fasta_dict = extract_low_coverage_sequences(fasta_reader, low_cov_regions)
    write_fasta_file(low_cov_fasta_dict, low_cov_fasta_out)

    # Step 5: Extract assembly without low-coverage regions and save to a dictionary
    logging.info("Extracting low-coverage regions from assembly.")
    fasta_reader = fasta_reader_yield(fasta_file_path)  # Resetting the fasta_reader generator
    splitted_fasta_dict = extract_splitted_fasta(fasta_reader, low_cov_regions, remove_low_coverage)
    write_fasta_file(splitted_fasta_dict, splitted_fasta_out)
    logging.info(f"Splitted assebmbly is written in {splitted_fasta_out}")
    
    # Step 6: Write tables
    write_tables(low_cov_regions, high_cov_regions, contig_stats_df, low_cov_file, high_cov_file, contig_stats_file)

    # Step 7: Write final statistics report
    logging.info("Final stats:")
    write_final_statistics_report(low_cov_regions, high_cov_regions, final_stats_file)


if __name__ == "__main__":
    main()
