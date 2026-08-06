process DECONTAM_STATS {
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy conda-forge::scipy"
    container "python:3.11-slim"

    publishDir "${params.outdir}/stats", mode: params.publish_dir_mode

    input:
    path long_table     // combined_abundance_long.tsv
    path metadata       // sample_metadata.tsv

    output:
    path "taxon_stats.tsv",       emit: table
    path "significant_taxa.tsv",  emit: significant
    path "sankey_flags.tsv",      emit: flags
    path "stats_summary.json",    emit: summary_json

    script:
    """
    decontam_stats.py \\
        --abundance ${long_table} \\
        --metadata ${metadata} \\
        --rank ${params.stats_rank} \\
        --min-reads ${params.stats_min_reads} \\
        --min-abundance ${params.stats_min_abundance} \\
        --alpha ${params.stats_alpha} \\
        --min-fold ${params.stats_min_fold} \\
        --bg-rate ${params.stats_bg_rate} \\
        --out-prefix .
    """

    stub:
    """
    printf "sample\\ttaxid\\tname\\treads\\trel_abundance\\tp_value\\tq_value\\tcall\\tmode\\n" > taxon_stats.tsv
    printf "SAMPLE01\\t562\\tEscherichia coli\\t300\\t0.30\\t0.001\\t0.004\\tdetected\\tno_control\\n" >> taxon_stats.tsv
    cp taxon_stats.tsv significant_taxa.tsv
    printf "sample\\ttaxid\\tcall\\n SAMPLE01\\t562\\tdetected\\n" > sankey_flags.tsv
    printf '{"mode":"no_control","n_significant":1}' > stats_summary.json
    """
}
