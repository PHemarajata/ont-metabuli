process MEGAHIT {
    tag "${meta.id}:${bin_label}"
    label 'process_assembly'

    conda "bioconda::megahit>=1.2.9"
    container "quay.io/biocontainers/megahit:1.2.9--haf24da9_8"

    publishDir "${params.outdir}/assembly/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), val(taxid), path(reads)

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}.assembly.fasta"), path(reads), emit: assembly
    path "${meta.id}.${bin_label}.megahit.log", emit: log, optional: true

    script:
    def out = "${meta.id}.${bin_label}.assembly.fasta"
    """
    n=\$(( \$(wc -l < ${reads[0]}) / 4 ))
    echo "[megahit] ${meta.id} ${bin_label}: \$n pairs (min ${params.min_bin_reads})" | tee ${meta.id}.${bin_label}.megahit.log

    if [ "\$n" -lt ${params.min_bin_reads} ]; then
        echo "[megahit] too few reads — skipping assembly" | tee -a ${meta.id}.${bin_label}.megahit.log
        : > ${out}
        exit 0
    fi

    # MEGAHIT is the memory-light choice for short-read metagenome assembly.
    # It refuses to run if the output dir exists, hence --out-dir on a fresh path.
    if megahit -1 ${reads[0]} -2 ${reads[1]} \\
            --min-contig-len ${params.min_contig_len} \\
            -t ${task.cpus} -o megahit_out >> ${meta.id}.${bin_label}.megahit.log 2>&1; then
        if [ -s megahit_out/final.contigs.fa ]; then
            cp megahit_out/final.contigs.fa ${out}
        else
            : > ${out}
        fi
    else
        echo "[megahit] assembly failed (no contigs produced)" | tee -a ${meta.id}.${bin_label}.megahit.log
        : > ${out}
    fi
    """

    stub:
    """
    printf '>k141_1\\nACGTACGTACGTACGTACGT\\n' > ${meta.id}.${bin_label}.assembly.fasta
    echo stub > ${meta.id}.${bin_label}.megahit.log
    """
}
