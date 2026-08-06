process POLISH {
    tag "${meta.id}:${bin_label}"
    label 'process_medium'

    conda "bioconda::racon=1.5.0 bioconda::minimap2=2.30"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), path(assembly), path(reads)

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}.polished.fasta"), emit: assembly

    script:
    def out = "${meta.id}.${bin_label}.polished.fasta"
    if (params.polisher == 'none')
        """
        cp ${assembly} ${out}
        """
    else if (params.polisher == 'medaka')
        // medaka lives in a separate env (python/samtools conflict); use it if set
        """
        if [ ! -s ${assembly} ]; then cp ${assembly} ${out}; exit 0; fi
        MED="${params.medaka_env ?: ''}"
        if [ -n "\$MED" ] && [ -x "\$MED/bin/medaka_consensus" ]; then
            "\$MED/bin/medaka_consensus" -i ${reads} -d ${assembly} -o medaka_out -t ${task.cpus} ${params.medaka_model ? "-m ${params.medaka_model}" : ''} \\
              && cp medaka_out/consensus.fasta ${out} || cp ${assembly} ${out}
        else
            echo "[polish] medaka env not found (--medaka_env); falling back to racon" >&2
            minimap2 -x map-ont -t ${task.cpus} ${assembly} ${reads} > ov.paf
            racon -t ${task.cpus} ${reads} ov.paf ${assembly} > ${out} || cp ${assembly} ${out}
        fi
        """
    else   // default: racon (fast, model-free)
        """
        if [ ! -s ${assembly} ]; then cp ${assembly} ${out}; exit 0; fi
        minimap2 -x map-ont -t ${task.cpus} ${assembly} ${reads} > ov.paf
        racon -t ${task.cpus} ${reads} ov.paf ${assembly} > ${out} || cp ${assembly} ${out}
        """

    stub:
    "cp ${assembly} ${meta.id}.${bin_label}.polished.fasta"
}
