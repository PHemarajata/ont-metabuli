process SUMMARY {
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas"
    container "python:3.11-slim"

    publishDir "${params.outdir}/summary", mode: params.publish_dir_mode

    input:
    path sankeys,   stageAs: 'sankey/*'
    path kronas,    stageAs: 'krona/*'
    path significant
    path summary_json
    path matrix
    path qc_jsons,  stageAs: 'qc/*'
    path host_stats, stageAs: 'host/*'

    output:
    path "index.html",        emit: html
    path "run_summary.tsv",   emit: table

    script:
    """
    make_summary.py \\
        --sankey-dir sankey \\
        --krona-dir krona \\
        --significant ${significant} \\
        --summary-json ${summary_json} \\
        --qc-dir qc \\
        --host-dir host \\
        --out index.html \\
        --table run_summary.tsv
    """

    stub:
    """
    echo "<html><body>summary</body></html>" > index.html
    printf "sample\\tmetric\\tvalue\\n" > run_summary.tsv
    """
}
