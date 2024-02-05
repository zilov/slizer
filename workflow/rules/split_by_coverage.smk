rule split_by_coverage:
    input:
        assembly = ASSEMBLY,
        coverage_bed = rules.get_coverage.output.coverage_bed
    conda: envs.bedtools
    threads: workflow.cores
    output:
        splitted_fasta = f"{OUTDIR}/{PREFIX}_contig_stats.tsv",
        split_stats = f"{OUTDIR}/{PREFIX}_final_statistics_report.txt",
        fasta_splitted = f"{OUTDIR}/{PREFIX}_splitted.fa",
    params:
        min_length = MIN_LENGTH,
        outdir = directory(OUTDIR),
        split_script = SPLIT_SCRIPT,
        prefix = PREFIX,
        zl = ZL_THRESHOLD,
        zh = ZH_THRESHOLD,
        tail = TAIL_THRESHOLD,
        mean_coverage_fraction = MEAN_COVERAGE_FRACTION,
        remove_low_coverage = REMOVE_LOW_COVERAGE
    shell: "{params.split_script} -a {input.assembly} -c {params.mean_coverage_fraction} -l {params.min_length} \
                -zl {params.zl} -zh {params.zh} -t {params.tail} -b {input.coverage_bed} -p {params.prefix} -o {params.outdir} -r {params.remove_low_coverage}" 