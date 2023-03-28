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