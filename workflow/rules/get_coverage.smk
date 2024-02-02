rule get_coverage:
    input:
        bam = BAM
    threads: workflow.cores
    conda: envs.align
    output:
        coverage_bed = f"{OUTDIR}/alignment/{PREFIX}_coverage.bed",
    params:
        outdir = directory(f"{OUTDIR}/alignment")
    shell: """
        bedtools genomecov -bga -split -ibam {input.bam} > {output.coverage_bed}
    """