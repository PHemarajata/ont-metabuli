process SPADES {
    tag "${meta.id}:${bin_label}"
    label 'process_assembly'

    conda "bioconda::spades=4.3.0"
    container "quay.io/biocontainers/spades:4.3.0--hde4eca7_1"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), val(taxid), path(reads)

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}.assembly.fasta"), path(reads), emit: assembly
    path "${meta.id}.${bin_label}.spades.log", emit: log, optional: true

    script:
    def out = "${meta.id}.${bin_label}.assembly.fasta"
    def mem = task.memory ? task.memory.toGiga() : 8
    """
    n=\$(( \$(wc -l < ${reads[0]}) / 4 ))
    echo "[spades] ${meta.id} ${bin_label}: \$n pairs (min ${params.min_bin_reads})" | tee ${meta.id}.${bin_label}.spades.log

    if [ "\$n" -lt ${params.min_bin_reads} ]; then
        echo "[spades] too few reads — skipping assembly" | tee -a ${meta.id}.${bin_label}.spades.log
        : > ${out}
        exit 0
    fi

    # metaSPAdes: higher contiguity than MEGAHIT at a substantially larger
    # memory/time cost. -m is a hard cap in GB so the task cannot exceed its
    # requested memory.
    if spades.py --meta -1 ${reads[0]} -2 ${reads[1]} \\
            -t ${task.cpus} -m ${mem} -o spades_out >> ${meta.id}.${bin_label}.spades.log 2>&1; then
        if [ -s spades_out/contigs.fasta ]; then
            # apply the same minimum contig length used for MEGAHIT
            awk -v m=${params.min_contig_len} '/^>/{keep=0; n=\$0; sub(/.*_length_/,"",n); sub(/_cov.*/,"",n); if (n+0>=m) keep=1} keep' \\
                spades_out/contigs.fasta > ${out} || cp spades_out/contigs.fasta ${out}
        else
            : > ${out}
        fi
    else
        echo "[spades] assembly failed (no contigs produced)" | tee -a ${meta.id}.${bin_label}.spades.log
        : > ${out}
    fi
    """

    stub:
    """
    printf '>NODE_1_length_500_cov_10\\nACGTACGTACGTACGTACGT\\n' > ${meta.id}.${bin_label}.assembly.fasta
    echo stub > ${meta.id}.${bin_label}.spades.log
    """
}
