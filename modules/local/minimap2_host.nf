process MINIMAP2_HOST {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::minimap2=2.30 bioconda::samtools=1.23"
    container "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0"

    publishDir "${params.outdir}/host_removed", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.endsWith('.flagstat') ? fn : null }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("${meta.id}*.hostremoved.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.host.flagstat"),         emit: flagstat

    script:
    def preset = meta.paired ? 'sr' : 'map-ont'
    if (meta.paired)
        """
        # Map the pair to the host, then keep only pairs where BOTH mates are
        # unmapped (-f 12). collate restores mate adjacency so samtools fastq
        # can emit a properly paired R1/R2 set.
        minimap2 -ax ${preset} -t ${task.cpus} -I ${params.minimap2_index_size} \\
            ${index} ${reads[0]} ${reads[1]} 2> ${meta.id}.minimap2.log \\
            | samtools view -@ ${task.cpus} -b -o ${meta.id}.tmp.bam -

        samtools flagstat -@ ${task.cpus} ${meta.id}.tmp.bam > ${meta.id}.host.flagstat

        samtools view -@ ${task.cpus} -b -f 12 -F 256 ${meta.id}.tmp.bam \\
            | samtools collate -@ ${task.cpus} -u -O - \\
            | samtools fastq -@ ${task.cpus} -n \\
                -1 ${meta.id}_R1.hostremoved.fastq.gz \\
                -2 ${meta.id}_R2.hostremoved.fastq.gz \\
                -0 /dev/null -s /dev/null -

        rm -f ${meta.id}.tmp.bam
        """
    else
        """
        # Map to host (ONT preset). Stream to a temp BAM, tally, then keep only
        # UNMAPPED reads (-f 4) as the decontaminated read set. Memory is bounded
        # by minimap2 -I; no full-length host BAM is retained.
        minimap2 -ax ${preset} -t ${task.cpus} -I ${params.minimap2_index_size} \\
            ${index} ${reads} 2> ${meta.id}.minimap2.log \\
            | samtools view -@ ${task.cpus} -b -o ${meta.id}.tmp.bam -

        samtools flagstat -@ ${task.cpus} ${meta.id}.tmp.bam > ${meta.id}.host.flagstat

        samtools fastq -@ ${task.cpus} -n -f 4 ${meta.id}.tmp.bam \\
            | gzip > ${meta.id}.hostremoved.fastq.gz

        rm -f ${meta.id}.tmp.bam
        """

    stub:
    if (meta.paired)
        """
        echo "" | gzip > ${meta.id}_R1.hostremoved.fastq.gz
        echo "" | gzip > ${meta.id}_R2.hostremoved.fastq.gz
        printf "100 + 0 in total\\n10 + 0 mapped (10.00%%)\\n" > ${meta.id}.host.flagstat
        """
    else
        """
        echo "" | gzip > ${meta.id}.hostremoved.fastq.gz
        printf "100 + 0 in total\\n10 + 0 mapped (10.00%%)\\n" > ${meta.id}.host.flagstat
        """
}
