rule align_paired:
    input:
        assembly = ASSEMBLY,
        fr = FR,
        rr = RR    
    threads: workflow.cores
    conda: envs.align
    output:
        alignment_bam = f"{OUTDIR}/alignment/{PREFIX}_sorted.bam",
        alignment_stats = f"{OUTDIR}/alignment/{PREFIX}_sorted.stats",
    params:
        outdir = directory(f"{OUTDIR}/alignment")
    shell: """
    minimap2 -ax sr {input.assembly} {input.fr} {input.rr} -t {threads} \
    | samtools view -b -u - \
    | samtools sort -@ {threads} -T tmp > {output.alignment_bam} \
    && samtools flagstat {output.alignment_bam} > {output.alignment_stats} \
    """