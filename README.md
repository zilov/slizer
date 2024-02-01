# slizer
Tool for cutting assembly contigs by low coverage regions, which could have misassembly/chimeric origins. Slizer is designed to:
- run and sort alignment of sequencing reads on draft assembly with [`minimap2`](https://github.com/lh3/minimap2) and [`samtools`](http://www.htslib.org/)
- then identify coverage regions with [`bedtools genomecov`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
- report low and high coverage regions in genome based on z-score thresholds
- split contigs by low-coverage regions. 

Z-score is calculated for each contig separately. By default z-score value for low coverage regions is -2.0, for high coverage 2.0. 
To define region as low or high coverage, value of coverage for each interval should be consistent with filters (by default):
- z-score lower than -2.0 or highier than 2.0 (`-zl` and `-zh` arguments)
- length of region (`-l` argument, > 500 bp by default)
- for low coverage regions, coverage is lower than mean genome coverage fraction (`-c` argument, 0.5 by default)
- region is not in tails of contig (`--tail_threshold`, 1000 by default)
  
After z-scores is calculated, assembly contigs is splitted by low coverage regions. Low coverage regions could be excluded from resulting splitted assembly with `-r` argument.

## Installation

To install the `slizer`, you need to have Python 3.6 or later installed. Additionally, you should install `snakemake` using `conda` or `mamba`:

```bash
conda/mamba install -c bioconda -c conda-forge -c defaults snakemake
```

During first run snakemake will install all dependencies for slizer, you do not need to set it manually.

## Usage
To use the slizer tool with HiFi reads, run the following command:

```bash
slizer.py -a <assembly_fasta> -hifi <hifi_reads.fq> -o <output_folder>
```
By default slizer uses 8 threads to run `minimap2`.

You could also use Nanopore, PacBio or pair-end Illumina reads to run slizer. 

Also, if you already have reads alignment file, you could generage `coverage.bed` file directly and provide it to slizer. Use `bedtools genomecov -bga -split -ibam {your_alignment.bam} > {output.coverage_bed}` than run slizer with `slizer -m coverage -a <assembly.fasta> -b <coverage.bam> -o ./slizer_results`.

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
* -c, --mean_coverage_fraction: Fraction of the global mean coverage to use as a threshold for filtering low intervals. Intervals with coverage below fraction will be used for further analysis. (default: 0.5)
* --tail_threshold: Size of the tail regions to exclude. (default: 1000)
* -zl, --z_threshold_low: Z-score threshold for low coverage. Should be negative float value.
* -zh, --z_threshold_high: Z-score threshold for high coverage.
* -r, --remove_low_coverage: removes low coverage regions from splitted assembly
* -o, --output_folder: Path to the output folder (default: ./slizer)
* -d, --debug, debug mode / dry run

## Example

```bash
slizer.py -m coverage -a assembly.fasta -b coverage.bed -l 100 -c 0.2 -zl -1.8 -zh 4.0 -r -o ./output -t 50
```

This command will process the input `assembly.fasta` and `coverage.bed` created by user to create output files in the ./output folder using 50 threads. The minimum length threshold for coverage intervals is set to 100, and the mean coverage fraction is set to 0.2, low-coverage threahold is -1.8, high coverage threshold is 4.0. Low coverage reguos will be removed from final splitted assembly file.

## Output Files

`slizer` generates five output files in the specified output folder (default is `./slizer`):

1. `alignment`: folder with minimap2 alignment and alignment stats files 
2. `*_regions.bed`: A BED file containing the extracted low coverage regions based on the specified length and coverage thresholds.
3. `*_masked.fasta`: A FASTA file containing the input assembly with low coverage regions replaced with N's.
4. `*_splitted.fasta`: A FASTA file containing the contigs from the masked assembly, split by N's.
5. `*.stats`: A text file containing statistics about the processed assembly, such as the number of regions removed, bases removed, and contigs split.
6. `*.log`: A log file containing information about the processing steps and any errors or warnings encountered during the run.

The output files are named using the assembly prefix,or prefix set by user.

Example:

```bash
alignment/
my_assembly_regions.bed
my_assembly_masked.fasta
my_assembly_splitted.fasta
my_assembly.stats
my_assembly.log
```
