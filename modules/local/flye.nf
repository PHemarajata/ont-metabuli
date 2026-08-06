process FLYE {
    tag "${meta.id}:${bin_label}"
    label 'process_assembly'

    conda "bioconda::flye=2.9.6"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), val(taxid), path(reads)

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}.assembly.fasta"), path(reads), emit: assembly
    path "${meta.id}.${bin_label}.flye.log",                                                       emit: log, optional: true

    script:
    def out = "${meta.id}.${bin_label}.assembly.fasta"
    """
    n=\$(( \$(wc -l < ${reads}) / 4 ))
    echo "[flye] ${meta.id} ${bin_label}: \$n reads (min ${params.min_bin_reads})" | tee ${meta.id}.${bin_label}.flye.log

    if [ "\$n" -lt ${params.min_bin_reads} ]; then
        echo "[flye] too few reads — skipping assembly" | tee -a ${meta.id}.${bin_label}.flye.log
        : > ${out}
        exit 0
    fi

    # metaFlye; --meta handles uneven coverage (single-taxon or diverse bins alike)
    if flye --${params.flye_read_type} ${reads} --meta \\
            --threads ${task.cpus} --out-dir flye_out >> ${meta.id}.${bin_label}.flye.log 2>&1; then
        if [ -s flye_out/assembly.fasta ]; then
            cp flye_out/assembly.fasta ${out}
            [ -f flye_out/assembly_info.txt ] && cp flye_out/assembly_info.txt ${meta.id}.${bin_label}.assembly_info.txt || true
        else
            : > ${out}
        fi
    else
        echo "[flye] assembly failed (no contigs produced)" | tee -a ${meta.id}.${bin_label}.flye.log
        : > ${out}
    fi
    """

    stub:
    """
    printf '>contig_1\\nACGTACGTACGTACGTACGT\\n' > ${meta.id}.${bin_label}.assembly.fasta
    echo stub > ${meta.id}.${bin_label}.flye.log
    """
}
