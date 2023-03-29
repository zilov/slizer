#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#@created: 28.03.2023
#@author: Danil Zilov
#@contact: zilov.d@gmail.com

import subprocess
import os
import argparse
import shutil
import logging
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
        logging.error(f"An error occurred while running bedtools maskfasta: {e}")
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

    with open(bed_file) as fh:
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

def validate_args(args: argparse.Namespace) -> None:
    if not os.path.isfile(args.assembly) or os.path.getsize(args.assembly) == 0:
        raise ValueError(f"The assembly file '{args.assembly}' does not exist or empty.")
    if not os.path.isfile(args.coverage_bed) or os.path.getsize(args.coverage_bed) == 0:
        raise ValueError(f"The coverage bed file '{args.coverage_bed}' does not exist or empty.")
    if args.min_length_threshold < 1:
        raise ValueError("The minimum length threshold should be a positive integer.")
    if args.max_coverage_threshold < 1:
        raise ValueError("The maximum coverage threshold should be a positive integer.")

def check_bedtools_availability():
    if not shutil.which("bedtools"):
        raise FileNotFoundError("The 'bedtools' executable was not found. Please ensure it is installed and available in your system PATH.")
    
def check_output_files(output_files: List[str]):
    for file in output_files:
        if not os.path.isfile(file):
            raise FileNotFoundError(f"An error occurred while generating the output file '{file}'. Please check the input files and parameters.")
        
def configure_logging(log_file):
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s]: %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s]: %(message)s"))

    # Add the file handler to the logger
    logging.getLogger().addHandler(file_handler)

def main(args: argparse.Namespace) -> Tuple[str, str, str, str]:
    
    # Convert input arguments to named tuples
    assembly = os.path.abspath(args.assembly)
    coverage_bed = os.path.abspath(args.coverage_bed)
    output_folder = os.path.abspath(args.output_folder)
    min_length_threshold = args.min_length_threshold
    max_coverage_threshold = args.max_coverage_threshold
    
    # Get the assembly prefix from the file name
    assembly_prefix = os.path.splitext(os.path.basename(assembly))[0]
    os.makedirs(output_folder, exist_ok=True)
    
    log_file = os.path.join(args.output_folder, f"{assembly_prefix}_l{args.min_length_threshold}_c{args.max_coverage_threshold}.log")
    configure_logging(log_file)
    
    logging.info(f"Length threshold is {min_length_threshold}, coverage threshold is {max_coverage_threshold}")
    logging.info(f"Assembly: {assembly}, coverage file: {coverage_bed}")
    

    # Set the output file paths with the assembly prefix and thresholds
    regions_bed = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_regions.bed")
    masked_fasta = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_masked.fasta")
    split_fasta = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}_splitted.fasta")
    run_stats = os.path.join(output_folder, f"{assembly_prefix}_l{min_length_threshold}_c{max_coverage_threshold}.stats")
    
    
    # Extract regions from the coverage BED file based on the coverage thresholds
    logging.info("Extracting low coverage regions by thresholds...")
    extract_coverage_gaps(min_length_threshold, max_coverage_threshold, coverage_bed, regions_bed)
    logging.info("Done!")
    
    # Mask the input assembly FASTA file with N's using the BED file
    try:
        logging.info("Replacing low coverage regions with N's...")
        insert_N_to_regions(regions_bed, assembly, masked_fasta)
        logging.info("Done!")
    except subprocess.CalledProcessError as e:
        logging.error(f"An error occurred while running bedtools maskfasta: {e}")
    
    # Split the masked FASTA file by N's and write the split sequences to a new FASTA file
    with open(split_fasta, "w") as out_fh:
        logging.info(f"Splitting {masked_fasta} by N's...")
        for header, seq in split_fasta_by_N(masked_fasta):
            out_fh.write(f">{header}\n{seq}\n")
        logging.info("Done!")
            
    # Write run stats
    with open(run_stats, 'w') as fw:
        fw.write(f"Length threshold is {min_length_threshold}, coverage threshold is {max_coverage_threshold}")
        fw.write(f"Assembly: {assembly}, coverage file: {coverage_bed}\n")
        bases, contigs, regions = count_bed_stats(regions_bed)
        fw.write(f"Stats:\nRegions removed:\t{regions}\n")
        fw.write(f"Bases removed:\t{bases}\n")
        fw.write(f"Contigs splitted:\t{contigs}\n")
        logging.info(f"Stats:")
        logging.info(f"Regions removed: {regions}")
        logging.info(f"Bases removed: {bases}")
        logging.info(f"Contigs splitted: {contigs}")
        
    # At the end of the main function, add the return statement
    
    output_files = [regions_bed, masked_fasta, split_fasta, run_stats]
    check_output_files(output_files)
    logging.info(f"All done! Output folder is: {output_folder}")
    
    return regions_bed, masked_fasta, split_fasta, run_stats


def run_script():
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Process assembly fasta and coverage bed files.')
    parser.add_argument('-a', '--assembly', type=str, help='Path to the assembly fasta file.', required=True)
    parser.add_argument('-b', '--coverage_bed', type=str, help='Path to the coverage bed file.')
    parser.add_argument('-l', '--min_length_threshold', type=int, default=100, help='Minimum length threshold for coverage gaps.')
    parser.add_argument('-c', '--max_coverage_threshold', type=int, default=10, help='Maximum coverage threshold for coverage gaps.')
    parser.add_argument('-o', '--output_folder', type=str, default="./asShredder", help='Path to the output folder.')
    
    args = parser.parse_args()
    validate_args(args)
    check_bedtools_availability()
    main(args)

if __name__ == "__main__":
    run_script()