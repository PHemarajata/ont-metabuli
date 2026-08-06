#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ont-metabuli
    Low-resource Oxford Nanopore metagenomics classification pipeline

      samplesheet ->  QC (nanoq)
                  ->  host removal (minimap2 map-ont, keep unmapped)
                  ->  classification (Metabuli, long-read mode)
                  ->  combine reports  ->  decontam / significance stats
                  ->  per-sample Sankey (+ Krona)  ->  HTML summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
//  Modules
// ---------------------------------------------------------------------------
include { CHOPPER          } from './modules/local/chopper.nf'
include { FASTP            } from './modules/local/fastp.nf'
include { MINIMAP2_INDEX   } from './modules/local/minimap2_index.nf'
include { MINIMAP2_HOST    } from './modules/local/minimap2_host.nf'
include { METABULI_CLASSIFY} from './modules/local/metabuli.nf'
include { COMBINE_REPORTS  } from './modules/local/combine_reports.nf'
include { DECONTAM_STATS   } from './modules/local/decontam_stats.nf'
include { SANKEY           } from './modules/local/sankey.nf'
include { SUMMARY          } from './modules/local/summary.nf'
include { ASSEMBLY         } from './subworkflows/local/assembly.nf'

// ---------------------------------------------------------------------------
//  Help
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    ================================================================
     ont-metabuli  v${workflow.manifest.version}
    ================================================================
    Usage:
      nextflow run . -profile standard \\
          --input samplesheet.csv \\
          --metabuli_db /path/to/metabuli_db \\
          --host_fasta /path/to/chm13v2.0.fa.gz \\
          --outdir results

    Samplesheet (CSV, header required):
      ONT (single file per sample):
        sample,fastq,sample_type,batch
        SAMPLE01,/abs/path/reads.fastq.gz,sample,run1
        NTC01,/abs/path/ntc.fastq.gz,negative_control,run1

      Illumina paired-end (add fastq_2, and pass --platform illumina):
        sample,fastq,fastq_2,sample_type,batch
        SAMPLE01,/abs/path/S1_R1.fastq.gz,/abs/path/S1_R2.fastq.gz,sample,run1

      * sample_type: 'sample' for biological samples; one of
        [${params.control_types}] for negative controls (optional).
      * batch: optional; used to pair controls to samples. Blank = one pool.

    Platform:
      --platform ont|illumina       default '${params.platform}'. Selects chopper vs fastp,
                                    minimap2 map-ont vs sr, Metabuli --seq-mode 3 vs 2,
                                    and Flye vs MEGAHIT/SPAdes. Set per RUN.

    Key params (see nextflow.config for all):
      --min_qual / --min_len        ONT QC thresholds (default ${params.min_qual} / ${params.min_len})
      --sr_min_qual / --sr_min_len  Illumina QC thresholds (default ${params.sr_min_qual} / ${params.sr_min_len})
      --metabuli_db                 Metabuli DB directory (required)
      --metabuli_max_ram            GiB cap for Metabuli (default ${params.metabuli_max_ram})
      --host_fasta                  host FASTA or .mmi (required unless --skip_host_removal)
      --stats_rank                  species | genus (default ${params.stats_rank})
      --skip_qc / --skip_host_removal
      --max_cpus / --max_memory / --max_time

    Optional assembly on binned reads (all OFF by default):
      --assemble_unclassified       de novo (Flye) on unclassified reads (novelty)
      --assemble_taxa '562,10239'   de novo per taxon; or 'significant' to use
                                    the taxa flagged by the stats step
      --reference_consensus 'T:/ref.fa;T2:/ref2.fa'
                                    map a taxon's reads to a reference -> consensus + variants
      --polisher racon|medaka|none  (default racon)   --classify_contigs
      --flye_read_type nano-hq|nano-raw|nano-corr     --min_bin_reads ${params.min_bin_reads}
    ================================================================
    """.stripIndent()
}

// ---------------------------------------------------------------------------
//  Samplesheet parsing
// ---------------------------------------------------------------------------
def parseSamplesheet(csvPath) {
    def controlSet = params.control_types.toString().toLowerCase().split(',').collect { it.trim() } as Set
    def illumina   = params.platform.toString().toLowerCase() == 'illumina'
    Channel
        .fromPath(csvPath, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) error "Samplesheet error: a row is missing 'sample'."
            if (!row.fastq)  error "Samplesheet error: sample '${row.sample}' is missing 'fastq'."
            def reads = [ file(row.fastq, checkIfExists: true) ]
            // 'fastq_2' turns a row into a read pair (Illumina)
            if (row.fastq_2 && row.fastq_2.toString().trim()) {
                if (!illumina)
                    error "Samplesheet error: sample '${row.sample}' has fastq_2 but --platform is '${params.platform}'. Use --platform illumina for paired-end data."
                reads << file(row.fastq_2, checkIfExists: true)
            } else if (illumina) {
                error "Samplesheet error: --platform illumina but sample '${row.sample}' has no 'fastq_2' column value. Paired-end reads are required; use --platform ont for single-file input."
            }
            def stype = (row.sample_type ?: 'sample').toString().trim().toLowerCase()
            def is_control = controlSet.contains(stype)
            def batch = (row.batch ?: 'all').toString().trim()
            def meta = [ id: row.sample.toString().trim(),
                         sample_type: is_control ? 'control' : 'sample',
                         raw_type: stype,
                         batch: batch ?: 'all',
                         is_control: is_control,
                         paired: reads.size() == 2 ]
            tuple(meta, reads)
        }
}

// ---------------------------------------------------------------------------
//  Workflow
// ---------------------------------------------------------------------------
workflow {

    if (params.help) { helpMessage(); return }

    // ---- validate required params ----
    if (!params.input)       error "Missing --input (samplesheet CSV). Run with --help."
    if (!params.metabuli_db) error "Missing --metabuli_db (Metabuli database directory)."
    if (!params.skip_host_removal && !params.host_fasta)
        error "Missing --host_fasta (host reference). Use --skip_host_removal to disable host filtering."

    def platform = params.platform.toString().toLowerCase()
    if (!(platform in ['ont', 'illumina']))
        error "--platform must be 'ont' or 'illumina' (got '${params.platform}')."
    is_illumina = platform == 'illumina'

    def qc_tool  = params.skip_qc ? 'SKIPPED' : (is_illumina ? 'fastp (paired-end)' : 'chopper (long read)')
    def asm_tool = (params.assembler ?: (is_illumina ? 'megahit' : 'flye'))

    // MEGAHIT's bioconda osx-arm64 build segfaults inside its 'assemble' step
    // (verified: exit -11 on any input), silently yielding empty assemblies.
    // The linux/amd64 container is fine, so this only affects native Apple Silicon.
    def on_mac_arm = System.getProperty('os.name')?.toLowerCase()?.contains('mac') &&
                     System.getProperty('os.arch') in ['aarch64', 'arm64']
    if (asm_tool == 'megahit' && on_mac_arm && !workflow.profile.contains('docker')) {
        log.warn "MEGAHIT's osx-arm64 build crashes during assembly and will produce empty contigs.\n" +
                 "         On Apple Silicon use '--assembler spades' (verified working) or run with '-profile docker'."
    }

    log.info """
    ================================================================
     ont-metabuli v${workflow.manifest.version}   |  profile: ${workflow.profile}
    ----------------------------------------------------------------
     platform     : ${platform}${is_illumina ? '  (paired-end)' : '  (long read)'}
     input        : ${params.input}
     QC           : ${qc_tool}
     metabuli_db  : ${params.metabuli_db}   (max-ram ${params.metabuli_max_ram} GiB, seq-mode ${params.metabuli_seq_mode ?: (is_illumina ? 2 : 3)})
     host removal : ${params.skip_host_removal ? 'SKIPPED' : "${params.host_fasta} (minimap2 ${is_illumina ? 'sr' : 'map-ont'})"}
     stats rank   : ${params.stats_rank}   (alpha ${params.stats_alpha})
     assembler    : ${asm_tool}${(params.assemble_unclassified || params.assemble_taxa || params.reference_consensus) ? '' : '  (assembly off)'}
     outdir       : ${params.outdir}
     resources    : ${params.max_cpus} cpus / ${params.max_memory} / ${params.max_time}
    ================================================================
    """.stripIndent()

    ch_reads = parseSamplesheet(params.input)

    // ---- 1. QC / read filtering (platform-specific) ----
    if (!params.skip_qc) {
        if (is_illumina) {
            FASTP(ch_reads)                 // paired-end trim + quality filter
            ch_filtered = FASTP.out.reads
            ch_qc_json  = FASTP.out.report
        } else {
            CHOPPER(ch_reads)               // long-read length + mean-quality filter
            ch_filtered = CHOPPER.out.reads
            ch_qc_json  = CHOPPER.out.report
        }
    } else {
        ch_filtered = ch_reads
        ch_qc_json  = Channel.empty()
    }

    // ---- 2. Host removal (minimap2, keep unmapped) ----
    if (!params.skip_host_removal) {
        ch_host = file(params.host_fasta, checkIfExists: true)
        if (params.host_fasta.toString().endsWith('.mmi')) {
            ch_index = Channel.value(ch_host)          // prebuilt index
        } else {
            MINIMAP2_INDEX(ch_host)
            ch_index = MINIMAP2_INDEX.out.index
        }
        MINIMAP2_HOST(ch_filtered, ch_index)
        ch_clean     = MINIMAP2_HOST.out.reads
        ch_host_stat = MINIMAP2_HOST.out.flagstat
    } else {
        ch_clean     = ch_filtered
        ch_host_stat = Channel.empty()
    }

    // ---- 3. Metabuli classification (long-read mode) ----
    ch_db = file(params.metabuli_db, checkIfExists: true, type: 'dir')
    METABULI_CLASSIFY(ch_clean, ch_db)

    // ---- 4. Combine per-sample reports into abundance tables ----
    ch_reports  = METABULI_CLASSIFY.out.report.map { meta, rep -> rep }.collect()
    ch_metadata = ch_reads
        .map { meta, fq -> "${meta.id}\t${meta.sample_type}\t${meta.batch}" }
        .collectFile(name: 'sample_metadata.tsv', newLine: true,
                     seed: 'sample\tsample_type\tbatch')
    COMBINE_REPORTS(ch_reports, ch_metadata)

    // ---- 5. Decontamination / significance statistics ----
    DECONTAM_STATS(COMBINE_REPORTS.out.long_table, ch_metadata)

    // ---- 6. Per-sample Sankey (coloured by stats flags when available) ----
    ch_flags = DECONTAM_STATS.out.flags   // single file: sample-taxon calls
    ch_for_sankey = METABULI_CLASSIFY.out.report.combine(ch_flags)
    SANKEY(ch_for_sankey)

    // ---- 6b. Optional assembly on binned reads (novelty / characterise / confirm) ----
    if (params.assemble_unclassified || params.assemble_taxa || params.reference_consensus) {
        // pair each sample's classified reads with its classification file
        ch_reads_class = ch_clean
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(METABULI_CLASSIFY.out.classifications.map { meta, cls -> tuple(meta.id, cls) })
            .map { id, meta, reads, cls -> tuple(meta, reads, cls) }
        ASSEMBLY(ch_reads_class, ch_db, DECONTAM_STATS.out.significant)
    }

    // ---- 7. HTML summary ----
    SUMMARY(
        SANKEY.out.html.map { meta, h -> h }.collect().ifEmpty([]),
        METABULI_CLASSIFY.out.krona.map { meta, k -> k }.collect().ifEmpty([]),
        DECONTAM_STATS.out.significant,
        DECONTAM_STATS.out.summary_json,
        COMBINE_REPORTS.out.matrix,
        ch_qc_json.map { meta, j -> j }.collect().ifEmpty([]),
        ch_host_stat.map { meta, f -> f }.collect().ifEmpty([])
    )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n✅ Pipeline completed. Open: ${params.outdir}/summary/index.html\n"
        : "\n❌ Pipeline failed. See .nextflow.log\n" )
}
