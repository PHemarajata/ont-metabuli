process ASSEMBLY_STATS {
    tag "${meta.id}:${bin_label}"
    label 'process_low'

    conda "bioconda::seqkit=2.13.0"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), path(assembly)

    output:
    path "${meta.id}.${bin_label}.assembly_stats.tsv", emit: stats

    script:
    """
    if [ -s ${assembly} ]; then
        seqkit stats -a -T ${assembly} \\
          | awk -v s='${meta.id}' -v b='${bin_label}' 'NR==1{print "sample\\tbin\\t"\$0} NR>1{print s"\\t"b"\\t"\$0}' \\
          > ${meta.id}.${bin_label}.assembly_stats.tsv
    else
        printf 'sample\\tbin\\tnote\\n%s\\t%s\\tno_assembly\\n' '${meta.id}' '${bin_label}' \\
          > ${meta.id}.${bin_label}.assembly_stats.tsv
    fi
    """

    stub:
    "printf 'sample\\tbin\\tfile\\nX\\t${bin_label}\\tstub\\n' > ${meta.id}.${bin_label}.assembly_stats.tsv"
}
