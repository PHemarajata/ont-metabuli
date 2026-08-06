process SANKEY {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::plotly"
    container "python:3.11-slim"

    publishDir "${params.outdir}/sankey", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(report), path(flags)

    output:
    tuple val(meta), path("${meta.id}.sankey.html"), emit: html

    script:
    """
    make_sankey.py \\
        --report ${report} \\
        --sample ${meta.id} \\
        --flags ${flags} \\
        --ranks ${params.sankey_ranks} \\
        --top-n ${params.sankey_top_n} \\
        --min-reads ${params.sankey_min_reads} \\
        --out ${meta.id}.sankey.html
    """

    stub:
    """
    echo "<html><body>Sankey ${meta.id}</body></html>" > ${meta.id}.sankey.html
    """
}
