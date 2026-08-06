process METABULI_EXTRACT {
    tag "${meta.id}:${bin_label}"
    label 'process_low'

    conda "bioconda::metabuli=1.2.0"
    container "quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0"

    input:
    tuple val(meta), path(reads), path(classifications), val(taxid), val(bin_label)
    path db

    output:
    tuple val(meta), val(bin_label), val(taxid), path("${meta.id}.${bin_label}*.fastq"), emit: reads

    script:
    // taxid -1 = unclassified reads (novelty mode)
    def seq_mode = meta.paired ? 2 : 3
    if (meta.paired)
        """
        mkdir -p ext
        # metabuli names outputs after the INPUT basenames, so stage canonical
        # names to make the results predictable: <sample>_R{1,2}_<taxid>.fq
        ln -sf ${reads[0]} q_R1.fastq.gz
        ln -sf ${reads[1]} q_R2.fastq.gz

        metabuli extract q_R1.fastq.gz q_R2.fastq.gz ${classifications} ${db} \\
            --tax-id ${taxid} --seq-mode ${seq_mode} --extract-format 2 --outdir ext || true

        # -1 becomes '-1' in the filename for the unclassified bin
        t=\$(echo "${taxid}" | sed 's/^-/-/')
        if [ -s "ext/q_R1_\${t}.fq" ] && [ -s "ext/q_R2_\${t}.fq" ]; then
            mv "ext/q_R1_\${t}.fq" ${meta.id}.${bin_label}_R1.fastq
            mv "ext/q_R2_\${t}.fq" ${meta.id}.${bin_label}_R2.fastq
        else
            # fall back to whatever was produced, in sorted (R1,R2) order
            set -- \$(ls ext/*.fq 2>/dev/null | sort)
            if [ -n "\${1:-}" ] && [ -n "\${2:-}" ]; then
                mv "\$1" ${meta.id}.${bin_label}_R1.fastq
                mv "\$2" ${meta.id}.${bin_label}_R2.fastq
            else
                : > ${meta.id}.${bin_label}_R1.fastq
                : > ${meta.id}.${bin_label}_R2.fastq
            fi
        fi
        echo "[extract] ${meta.id} ${bin_label} (taxid ${taxid}): \$(( \$(wc -l < ${meta.id}.${bin_label}_R1.fastq) / 4 )) pairs" >&2
        """
    else
        """
        mkdir -p ext
        metabuli extract ${reads} ${classifications} ${db} \\
            --tax-id ${taxid} --seq-mode ${seq_mode} --extract-format 2 --outdir ext || true

        found=\$(ls ext/*.fq 2>/dev/null | head -1 || true)
        if [ -n "\$found" ]; then
            mv "\$found" ${meta.id}.${bin_label}.fastq
        else
            : > ${meta.id}.${bin_label}.fastq   # empty bin -> downstream skips
        fi
        echo "[extract] ${meta.id} ${bin_label} (taxid ${taxid}): \$(( \$(wc -l < ${meta.id}.${bin_label}.fastq) / 4 )) reads" >&2
        """

    stub:
    if (meta.paired)
        """
        printf '@r/1\\nACGT\\n+\\n!!!!\\n' > ${meta.id}.${bin_label}_R1.fastq
        printf '@r/2\\nACGT\\n+\\n!!!!\\n' > ${meta.id}.${bin_label}_R2.fastq
        """
    else
        """
        printf '@r\\nACGT\\n+\\n!!!!\\n' > ${meta.id}.${bin_label}.fastq
        """
}
