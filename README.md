# asShredder
Tool for cutting assembly contigs by low coverage regions, which could have chimeric origins. 

## Installation

To install the `asShredder` tool, you need to have Python 3.6 or later installed. Additionally, you should install [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) using `conda` or `mamba`:

```bash
conda/mamba install -c bioconda bedtools
```

## Usage
To use the asShredder tool, run the following command:

```bash
python asShredder.py -a <assembly_fasta> -b <coverage_bed> [-l <min_length_threshold>] [-c <max_coverage_threshold>] [-o <output_folder>]
```

**Parameters:**

* -a, --assembly: Path to the assembly fasta file (required)
* -b, --coverage_bed: Path to the coverage bed file (required)
* -l, --min_length_threshold: Minimum length threshold for coverage gaps (default: 100)
* -c, --max_coverage_threshold: Maximum coverage threshold for coverage gaps (default: 10)
* -o, --output_folder: Path to the output folder (default: ./asShredder)

## Example

```bash
python asShredder.py -a assembly.fasta -b coverage.bed -l 100 -c 10 -o ./output
```

This command will process the input assembly.fasta and coverage.bed files and create output files in the ./output folder. The minimum length threshold for coverage gaps is set to 100, and the maximum coverage threshold is set to 10.

## Output Files

`asShredder` generates five output files in the specified output folder (default is `./asShredder`):

1. `*_regions.bed`: A BED file containing the extracted low coverage regions based on the specified length and coverage thresholds.
2. `*_masked.fasta`: A FASTA file containing the input assembly with low coverage regions replaced with N's.
3. `*_splitted.fasta`: A FASTA file containing the contigs from the masked assembly, split by N's.
4. `*.stats`: A text file containing statistics about the processed assembly, such as the number of regions removed, bases removed, and contigs split.
5. `*.log`: A log file containing information about the processing steps and any errors or warnings encountered during the run.

The output files are named using the assembly prefix, followed by the specified length and coverage thresholds.

Example:

```bash
my_assembly_l100_c10_regions.bed
my_assembly_l100_c10_masked.fasta
my_assembly_l100_c10_splitted.fasta
my_assembly_l100_c10.stats
my_assembly_l100_c10.log
```