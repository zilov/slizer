# asShredder
Tool for cutting assembly contigs by low coverage regions, which could have chimeric origins. asShredder is designed to run and sort alignment of sequencing reads on draft assembly with [`minimap2`](https://github.com/lh3/minimap2) and [`samtools`](http://www.htslib.org/), then identify coverage regions with [`bedtools genomecov`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html) and split contigs by length/coverage thresholds.

## Installation

To install the `asShredder` tool, you need to have Python 3.6 or later installed. Additionally, you should install `snakemake` using `conda` or `mamba`:

```bash
conda/mamba install -c bioconda -c conda-forge -c defaults snakemake
```

During first run snakemake will install all dependencies for asShredder, you do not need to set it manually.

## Usage
To use the asShredder tool with HiFi reads, run the following command:

```bash
asShredder.py -a <assembly_fasta> -hifi <hifi_reads.fq> -o <output_folder>
```
By default asShredder uses 500bp length threshold, 3 coverage thresholds and 8 threads to run `minimap2`.

You could also use Nanopore, PacBio or pair-end Illumina reads to run asShredder. 

Also, if you already have reads alignment file, you could generage `coverage.bed` file directly and provide it to asShredder. Use `bedtools genomecov -bga -split -ibam {your_alignment.bam} > {output.coverage_bed}` than run asShredder with `asShredder -m coverage -a <assembly.fasta> -b <coverage.bam> -o ./asshredder_results`.

**Parameters:**

* -m, --mode, mode to use {hifi,paired_reads,nanopore,pacbio,coverage} (hifi is default)
* -a, --assembly: Path to the assembly fasta file (required)
* -hifi, Path to HiFi reads
* -nanopore, Path to HiFi reads
* -pacbio, Path to PacBio reads
* -1, --forward_reads, Path to the forward reads fastq file
* -2, --reverse_reads, Path to the reverse reads fastq file
* -b, --coverage_bed: Path to the coverage bed file
* -p, --prefix,output files prefix (assembly file prefix by default)
* -l, --min_length_threshold: Minimum length threshold for coverage gaps (default: 500)
* -c, --max_coverage_threshold: Maximum coverage threshold for coverage gaps (default: 3)
* -o, --output_folder: Path to the output folder (default: ./asShredder)
* -d, --debug, debug mode / dry run

## Example

```bash
asShredder.py -m coverage -a assembly.fasta -b coverage.bed -l 100 -c 10 -o ./output -t 50
```

This command will process the input `assembly.fasta` and `coverage.bed` created by user to create output files in the ./output folder using 50 threads. The minimum length threshold for coverage gaps is set to 100, and the maximum coverage threshold is set to 10.

## Output Files

`asShredder` generates five output files in the specified output folder (default is `./asShredder`):

1. `alignment`: folder with minimap2 alignment and alignment stats files 
2. `*_regions.bed`: A BED file containing the extracted low coverage regions based on the specified length and coverage thresholds.
3. `*_masked.fasta`: A FASTA file containing the input assembly with low coverage regions replaced with N's.
4. `*_splitted.fasta`: A FASTA file containing the contigs from the masked assembly, split by N's.
5. `*.stats`: A text file containing statistics about the processed assembly, such as the number of regions removed, bases removed, and contigs split.
6. `*.log`: A log file containing information about the processing steps and any errors or warnings encountered during the run.

The output files are named using the assembly prefix, followed by the specified length and coverage thresholds.

Example:

```bash
alignment/
my_assembly_l100_c10_regions.bed
my_assembly_l100_c10_masked.fasta
my_assembly_l100_c10_splitted.fasta
my_assembly_l100_c10.stats
my_assembly_l100_c10.log
```