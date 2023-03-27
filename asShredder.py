#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#@created: 28.03.2023
#@author: Danil Zilov
#@contact: zilov.d@gmail.com

import subprocess
import os
import argparse
from typing import List, Tuple

def extract_coverage_gaps(min_length_threshold, max_coverage_threshold, coverage_bed, output_bed):
    command = "less %s | awk '{if ($4 < %s && ($3 - $2) > %s) print $0}' > %s" % (coverage_bed, max_coverage_threshold, min_length_threshold, output_bed)
    subprocess.run(command, shell=True)
    
def insert_N_to_regions(regions_bed, input_fasta, output_fasta):    
    # Run the bedtools maskfasta command using subprocess.run
    try:
        command = ["bedtools", "maskfasta", "-fi", input_fasta, "-fo", output_fasta, "-bed", regions_bed]
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bedtools maskfasta: {e}")
        raise e

def split_fasta_by_N(path_to_fasta_file):
    header = None
    with open(path_to_fasta_file) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                index = 0   
                if header:
                    # Split the sequence by N's
                    
                    seq_splits = filter(None, "".join(seq).split('N'))
                    # Yield the split sequences with new headers
                    for j, seq_split in enumerate(seq_splits):
                        seq_header = f"{header}_{index}"
                        yield seq_header, seq_split
                        index += 1
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            # Split the sequence by N's
            seq_splits = filter(None, "".join(seq).split('N'))
            # Yield the split sequences with new headers
            for j, seq_split in enumerate(seq_splits):
                seq_header = f"{header}_{index}"
                yield seq_header, seq_split
                index += 1

                
def count_bed_stats(bed_file):
    num_bases_deleted = 0
    num_contigs_split = 0
    num_regions_deleted = 0
    
    contigs_split = set()

    with open(bed_file, "r") as fh:
        for line in fh:
            fields = line.strip().split("\t")
            contig = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            if contig not in contigs_split:
                contigs_split.add(contig)
            num_bases_deleted += (end - start)
            num_regions_deleted += 1
        num_contigs_split = len(contigs_split)
                
    return num_bases_deleted, num_contigs_split, num_regions_deleted

def main(args: argparse.Namespace) -> Tuple[str, str, str, str]:
    
    # Convert input arguments to named tuples
    assembly = args.assembly
    coverage_bed = args.coverage_bed
    output_folder = args.output_folder
    min_length_threshold = args.min_length_threshold
    max_coverage_threshold = args.max_coverage_threshold

    # Get the assembly prefix from the file name
    assembly_prefix = os.path.splitext(os.path.basename(assembly))[0]
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Set the output file paths with the assembly prefix and thresholds
    regions_bed = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_regions.bed")
    masked_fasta = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_masked.fasta")
    split_fasta = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_splitted.fasta")
    run_stats = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}.stats")
    
    
    # Extract regions from the coverage BED file based on the coverage thresholds
    extract_coverage_gaps(min_length_threshold, max_coverage_threshold, coverage_bed, regions_bed)
    
    # Mask the input assembly FASTA file with N's using the BED file
    try:
        insert_N_to_regions(regions_bed, assembly, masked_fasta)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bedtools maskfasta: {e}")
    
    # Split the masked FASTA file by N's and write the split sequences to a new FASTA file
    with open(split_fasta, "w") as out_fh:
        for header, seq in split_fasta_by_N(masked_fasta):
            out_fh.write(f">{header}\n{seq}\n")
            
    # Write run stats
    with open(run_stats, 'w') as fw:
        fw.write(f"Length threshold is {min_length_threshold}, coverage threshold is {max_coverage_threshold}\n\n")
        bases, contigs, regions = count_bed_stats(regions_bed)
        fw.write(f"Stats:\nRegions removed: {regions}\n")
        fw.write(f"Bases removed: {bases}\n")
        fw.write(f"Contigs splitted: {contigs}\n")
    # At the end of the main function, add the return statement
    return regions_bed, masked_fasta, split_fasta, run_stats

if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Process assembly fasta and coverage bed files.')
    parser.add_argument('-a', '--assembly', type=str, help='Path to the assembly fasta file.', required=True)
    parser.add_argument('-b', '--coverage_bed', type=str, help='Path to the coverage bed file.')
    parser.add_argument('-l', '--min_length_threshold', type=int, help='Minimum length threshold for coverage gaps.')
    parser.add_argument('-c', '--max_coverage_threshold', type=int, help='Maximum coverage threshold for coverage gaps.')
    parser.add_argument('-o', '--output_folder', type=str, help='Path to the output folder.')
    
    args = parser.parse_args()
    main(args)