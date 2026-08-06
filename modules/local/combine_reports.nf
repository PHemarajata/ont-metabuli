process COMBINE_REPORTS {
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy"
    container "python:3.11-slim"

    publishDir "${params.outdir}/abundance", mode: params.publish_dir_mode

    input:
    path reports        // all *_report.tsv
    path metadata       // sample_metadata.tsv

    output:
    path "combined_abundance_long.tsv",       emit: long_table
    path "abundance_matrix_*.tsv",            emit: matrix
    path "read_accounting.tsv",               emit: accounting

    script:
    """
    combine_reports.py \\
        --reports ${reports} \\
        --metadata ${metadata} \\
        --rank ${params.stats_rank} \\
        --out-prefix combined
    """

    stub:
    """
    printf "sample\\ttaxid\\trank\\tname\\ttaxon_reads\\tclade_reads\\ttotal_classified\\trel_abundance\\n" > combined_abundance_long.tsv
    printf "SAMPLE01\\t562\\tS\\tEscherichia coli\\t300\\t300\\t1000\\t0.30\\n" >> combined_abundance_long.tsv
    printf "taxid\\tname\\tSAMPLE01\\n562\\tEscherichia coli\\t300\\n" > abundance_matrix_species.tsv
    printf "sample\\ttotal_classified\\n" > read_accounting.tsv
    """
}
