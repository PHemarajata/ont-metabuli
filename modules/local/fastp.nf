process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp>=0.23"
    container "quay.io/biocontainers/fastp:1.3.6--h43da1c4_0"

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode,
        saveAs: { fn -> (fn.endsWith('.json') || fn.endsWith('.html')) ? fn : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.filtered.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.qc.json"),                  emit: report
    path  "${meta.id}.fastp.html",                                emit: html

    script:
    def (r1, r2) = [ reads[0], reads[1] ]
    def adapters = params.adapter_fasta ? "--adapter_fasta ${file(params.adapter_fasta)}" : '--detect_adapter_for_pe'
    """
    # Paired-end trimming + quality filtering. fastp auto-detects adapters for
    # PE data unless an adapter FASTA is supplied.
    fastp \\
        --in1 ${r1} --in2 ${r2} \\
        --out1 ${meta.id}_R1.filtered.fastq.gz \\
        --out2 ${meta.id}_R2.filtered.fastq.gz \\
        ${adapters} \\
        --qualified_quality_phred ${params.sr_min_qual} \\
        --length_required ${params.sr_min_len} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        ${params.fastp_extra} \\
        2> ${meta.id}.fastp.log

    # Normalise fastp's JSON into the small schema the summary step reads.
    # Done with awk (no python/jq assumption about the container): find the
    # section, then take the first "total_reads" line inside it.
    get_total () {
        awk -v sec="\\"\$1\\"" 'index(\$0, sec){f=1} f && /"total_reads"/ {gsub(/[^0-9]/,""); print; exit}' \\
            ${meta.id}.fastp.json
    }
    before=\$(get_total before_filtering); before=\${before:-0}
    after=\$(get_total after_filtering);   after=\${after:-0}

    printf '{"sample":"%s","input_reads":%s,"reads":%s,"min_qual":%s,"min_len":%s}\\n' \\
        "${meta.id}" "\$before" "\$after" "${params.sr_min_qual}" "${params.sr_min_len}" \\
        > ${meta.id}.qc.json
    """

    stub:
    """
    echo "" | gzip > ${meta.id}_R1.filtered.fastq.gz
    echo "" | gzip > ${meta.id}_R2.filtered.fastq.gz
    printf '{"sample":"%s","input_reads":100,"reads":90}' "${meta.id}" > ${meta.id}.qc.json
    echo "<html>fastp</html>" > ${meta.id}.fastp.html
    """
}
