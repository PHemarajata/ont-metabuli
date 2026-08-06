process MINIMAP2_INDEX {
    tag "${host_fasta.name}"
    label 'process_medium'

    conda "bioconda::minimap2=2.30"
    container "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"

    publishDir "${params.outdir}/host_index", mode: params.publish_dir_mode

    input:
    path host_fasta

    output:
    path "host.*.mmi", emit: index

    script:
    // the index must be built with the same preset used for mapping
    def preset = params.platform.toString().toLowerCase() == 'illumina' ? 'sr' : 'map-ont'
    """
    # Preset-matched index (${preset}). Built once, reused across samples.
    minimap2 -x ${preset} -d host.${preset}.mmi -I ${params.minimap2_index_size} \\
        -t ${task.cpus} ${host_fasta}
    """

    stub:
    """
    touch host.${params.platform.toString().toLowerCase() == 'illumina' ? 'sr' : 'map-ont'}.mmi
    """
}
