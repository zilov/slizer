rule get_coverage_pacbio:
    input:
        assembly = ASSEMBLY,
        pacbio = PACBIO,    
    threads: workflow.cores
    conda: envs.align
    output:
        alignment_bam = f"{OUTDIR}/alignment/{PREFIX}_sorted.bam",
        alignment_stats = f"{OUTDIR}/alignment/{PREFIX}_sorted.stats",
        coverage_bed = f"{OUTDIR}/alignment/{PREFIX}_coverage.bed",
    params:
        outdir = directory(f"{OUTDIR}/alignment")
    shell: """
    minimap2 -ax map-pb {input.assembly} {input.pacbio} -t {threads} \
    | samtools view -b -u - \
    | samtools sort -@ {threads} -T tmp > {output.alignment_bam} \
    && samtools flagstat {output.alignment_bam} > {output.alignment_stats} \
    && bedtools genomecov -bga -split -ibam {output.alignment_bam} > {output.coverage_bed}
    """