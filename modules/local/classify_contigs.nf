process CLASSIFY_CONTIGS {
    tag "${meta.id}:${bin_label}"
    label 'process_metabuli'

    conda "bioconda::metabuli=1.2.0"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), path(assembly)
    path db

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}_contigs_report.tsv"), emit: report, optional: true

    when:
    // re-scans the (large) DB; only worth it to name novel/unclassified contigs
    params.classify_contigs

    script:
    """
    if [ ! -s ${assembly} ]; then
        echo "[classify_contigs] empty assembly — skipping" >&2
        exit 0
    fi
    mkdir -p out
    metabuli classify --seq-mode 3 --threads ${task.cpus} --max-ram ${params.metabuli_max_ram} \\
        ${assembly} ${db} out ${meta.id}_${bin_label}_contigs
    cp out/${meta.id}_${bin_label}_contigs_report.tsv ${meta.id}.${bin_label}_contigs_report.tsv
    """

    stub:
    """
    printf '100.0\\t1\\t1\\tspecies\\t562\\tEscherichia coli\\n' > ${meta.id}.${bin_label}_contigs_report.tsv
    """
}
