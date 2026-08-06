process REFERENCE_CONSENSUS {
    tag "${meta.id}:${bin_label}"
    label 'process_medium'

    conda "bioconda::minimap2=2.30 bioconda::samtools=1.23 bioconda::bcftools=1.24"

    publishDir "${params.outdir}/reference_consensus/${meta.id}/${bin_label}", mode: params.publish_dir_mode

    input:
    tuple val(meta), val(bin_label), path(reads), path(reference)

    output:
    tuple val(meta), val(bin_label), path("${meta.id}.${bin_label}.consensus.fasta"), emit: consensus, optional: true
    path "${meta.id}.${bin_label}.variants.vcf.gz",                                    emit: variants,  optional: true
    path "${meta.id}.${bin_label}.coverage.txt",                                       emit: coverage,  optional: true

    script:
    def pfx    = "${meta.id}.${bin_label}"
    def preset = meta.paired ? 'sr' : 'map-ont'
    def query  = meta.paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    def first  = meta.paired ? "${reads[0]}" : "${reads}"
    """
    if [ ! -s ${first} ]; then echo "[refcons] no reads for ${bin_label}" >&2; exit 0; fi

    # bcftools/samtools need a plain (faidx-able) reference
    REF=${reference}
    case "${reference}" in *.gz) gunzip -c ${reference} > ref.fa; REF=ref.fa;; esac
    samtools faidx "\$REF"

    minimap2 -ax ${preset} -t ${task.cpus} "\$REF" ${query} \\
        | samtools sort -@ ${task.cpus} -o aln.bam
    samtools index aln.bam

    samtools coverage aln.bam > ${pfx}.coverage.txt

    # simple consensus from the pileup
    samtools consensus -f fasta -o ${pfx}.consensus.fasta aln.bam || true

    # variants vs the reference
    bcftools mpileup -f "\$REF" aln.bam -Ou 2>/dev/null \\
        | bcftools call -mv -Oz -o ${pfx}.variants.vcf.gz 2>/dev/null || true
    [ -f ${pfx}.variants.vcf.gz ] && bcftools index ${pfx}.variants.vcf.gz || true

    echo "[refcons] ${meta.id} ${bin_label}: consensus + variants vs \$(basename ${reference})" >&2
    """

    stub:
    """
    printf '>consensus\\nACGT\\n' > ${meta.id}.${bin_label}.consensus.fasta
    echo "" | gzip > ${meta.id}.${bin_label}.variants.vcf.gz
    printf '#rname\\tstartpos\\tendpos\\tcoverage\\n' > ${meta.id}.${bin_label}.coverage.txt
    """
}
