process CHOPPER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::chopper=0.13.0"
    container "quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0"

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.json') ? fn : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.qc.json"),           emit: report

    script:
    """
    # chopper: length + mean-quality filter for ONT reads (Rust, streaming).
    chopper \\
        --quality ${params.min_qual} \\
        --minlength ${params.min_len} \\
        --threads ${task.cpus} \\
        --input ${reads} \\
        2> ${meta.id}.chopper.log \\
        | gzip > ${meta.id}.filtered.fastq.gz

    # distil the "Kept X reads out of Y reads" line into a small JSON report
    kept=\$(awk '/Kept/{print \$2}'  ${meta.id}.chopper.log); kept=\${kept:-0}
    total=\$(awk '/Kept/{print \$6}' ${meta.id}.chopper.log); total=\${total:-0}
    printf '{"sample":"%s","input_reads":%s,"reads":%s,"min_qual":%s,"min_len":%s}\\n' \\
        "${meta.id}" "\$total" "\$kept" "${params.min_qual}" "${params.min_len}" \\
        > ${meta.id}.qc.json
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.filtered.fastq.gz
    printf '{"sample":"%s","input_reads":100,"reads":90}' "${meta.id}" > ${meta.id}.qc.json
    """
}
