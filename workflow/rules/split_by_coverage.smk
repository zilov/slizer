rule split_by_coverage:
    input:
        assembly = ASSEMBLY,
        coverage_bed = COVERAGE_BED
    threads: workflow.cores
    conda: envs.align
    output:
        splitted_fasta = f"{OUTDIR}/{PREFIX}_l{MIN_LENGTH}_c{MAX_COVERAGE}_splitted.fasta",
        split_stats = f"{OUTDIR}/{PREFIX}_l{MIN_LENGTH}_c{MAX_COVERAGE}.stats",
    params:
        min_length = MIN_LENGTH,
        max_coverage = MAX_COVERAGE,
        outdir = OUTDIR,
        split_script = SPLIT_SCRIPT,
        prefix = PREFIX
    shell: "{params.split_script} \
    -a {input.assembly} -b {input.coverage_bed} \
    -l {params.min_length} -c {params.max_coverage} \
    -p {params.prefix} -o {params.outdir} "