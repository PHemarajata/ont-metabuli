process METABULI_CLASSIFY {
    tag "$meta.id"
    label 'process_metabuli'

    // On Apple Silicon prefer the native universal binary on PATH ('standard'
    // profile). conda/container are provided for x86 hosts / portability.
    conda "bioconda::metabuli=1.2.0"
    container "quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0"

    publishDir "${params.outdir}/metabuli/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path("${meta.id}_report.tsv"),          emit: report
    tuple val(meta), path("${meta.id}_classifications.tsv"), emit: classifications, optional: true
    tuple val(meta), path("${meta.id}_krona.html"),          emit: krona,           optional: true

    script:
    // seq-mode: 3 = long read, 2 = paired-end, 1 = single-end short read
    def seq_mode  = params.metabuli_seq_mode ?: (meta.paired ? 2 : 3)
    def min_score = params.metabuli_min_score > 0 ? "--min-score ${params.metabuli_min_score}" : ''
    // --precise 1 is Metabuli's short-read preset
    def precise   = (meta.paired && params.metabuli_precise) ? '--precise 1' : ''
    def query     = meta.paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    """
    mkdir -p out
    metabuli classify \\
        --seq-mode ${seq_mode} \\
        --threads ${task.cpus} \\
        --max-ram ${params.metabuli_max_ram} \\
        ${precise} ${min_score} ${params.metabuli_extra} \\
        ${query} ${db} out ${meta.id}

    # surface the three canonical outputs at the task root
    cp out/${meta.id}_report.tsv          ${meta.id}_report.tsv
    [ -f out/${meta.id}_classifications.tsv ] && cp out/${meta.id}_classifications.tsv ${meta.id}_classifications.tsv || true
    [ -f out/${meta.id}_krona.html ]          && cp out/${meta.id}_krona.html          ${meta.id}_krona.html          || true
    """

    stub:
    """
    cat > ${meta.id}_report.tsv <<'EOF'
100.0\t1000\t0\tR\t1\troot
95.0\t950\t0\tD\t2\t  Bacteria
40.0\t400\t50\tP\t1224\t    Proteobacteria
35.0\t350\t20\tG\t561\t      Escherichia
30.0\t300\t300\tS\t562\t        Escherichia coli
5.0\t50\t50\tS\t1280\t      Staphylococcus aureus
5.0\t50\t0\tU\t0\t  unclassified
EOF
    echo "<html>krona ${meta.id}</html>" > ${meta.id}_krona.html
    echo "classified ${meta.id}" > ${meta.id}_classifications.tsv
    """
}
