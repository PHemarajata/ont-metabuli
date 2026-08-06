//
//  Optional assembly on Metabuli-binned reads. Three modes, any combination:
//    A  --assemble_unclassified          de novo on unclassified reads (novelty)
//    B  --assemble_taxa '562,10239'|significant   de novo per taxon (characterise)
//    C  --reference_consensus 'taxid:/ref.fa;...' reference consensus + variants
//
//  All bins come from `metabuli extract`; a single EXTRACT step feeds both the
//  de novo path (Flye -> polish -> stats [-> classify contigs]) and the
//  reference-consensus path (minimap2 -> samtools consensus + bcftools).
//

include { METABULI_EXTRACT    } from '../../modules/local/metabuli_extract.nf'
include { FLYE                } from '../../modules/local/flye.nf'
include { MEGAHIT             } from '../../modules/local/megahit.nf'
include { SPADES              } from '../../modules/local/spades.nf'
include { POLISH              } from '../../modules/local/polish.nf'
include { ASSEMBLY_STATS      } from '../../modules/local/assembly_stats.nf'
include { CLASSIFY_CONTIGS    } from '../../modules/local/classify_contigs.nf'
include { REFERENCE_CONSENSUS } from '../../modules/local/reference_consensus.nf'

workflow ASSEMBLY {
    take:
    ch_reads_class    // tuple(meta, reads, classifications)
    ch_db             // path: metabuli DB dir
    ch_significant    // path: significant_taxa.tsv (for --assemble_taxa significant)

    main:
    def extract_inputs = Channel.empty()

    // ---- Mode A: unclassified reads (tax-id -1) ----
    if (params.assemble_unclassified) {
        extract_inputs = extract_inputs.mix(
            ch_reads_class.map { meta, reads, cls -> tuple(meta, reads, cls, -1, 'unclassified') }
        )
    }

    // ---- Mode B: per-taxon (explicit list or 'significant') ----
    if (params.assemble_taxa) {
        if (params.assemble_taxa.toString() == 'significant') {
            def sig = ch_significant
                .splitCsv(header: true, sep: '\t')
                .map { row -> tuple(row.sample.toString(), row.taxid.toString()) }
                .unique()
            def byid = ch_reads_class.map { meta, reads, cls -> tuple(meta.id, meta, reads, cls) }
            extract_inputs = extract_inputs.mix(
                sig.combine(byid, by: 0)
                   .map { id, taxid, meta, reads, cls -> tuple(meta, reads, cls, taxid, "taxon_${taxid}") }
            )
        } else {
            def taxids = params.assemble_taxa.toString().split(',').collect { it.trim() }.findAll { it }
            extract_inputs = extract_inputs.mix(
                ch_reads_class.flatMap { meta, reads, cls ->
                    taxids.collect { t -> tuple(meta, reads, cls, t, "taxon_${t}") } }
            )
        }
    }

    // ---- Mode C: reference consensus (extract taxon, then map to ref) ----
    def ch_refs = Channel.empty()
    if (params.reference_consensus) {
        def pairs = params.reference_consensus.toString().split(';')
            .collect { it.trim() }.findAll { it }
            .collect { p ->
                def i = p.indexOf(':')
                if (i < 0) error "Bad --reference_consensus entry '${p}' (want 'taxid:/path/ref.fa')"
                tuple(p.substring(0, i).trim(), file(p.substring(i + 1).trim(), checkIfExists: true))
            }
        ch_refs = Channel.fromList(pairs)                       // (taxid, reffile)
        extract_inputs = extract_inputs.mix(
            ch_reads_class.combine(Channel.fromList(pairs.collect { it[0] }))
                .map { meta, reads, cls, taxid -> tuple(meta, reads, cls, taxid, "ref_${taxid}") }
        )
    }

    // ---- single extraction step for every requested bin ----
    METABULI_EXTRACT(extract_inputs, ch_db)

    // route: reference-consensus bins vs de novo assembly bins
    def routed = METABULI_EXTRACT.out.reads.branch { meta, bin_label, taxid, reads ->
        refcons: bin_label.toString().startsWith('ref_')
        assemble: true
    }

    // ---- de novo path: assembler chosen by platform (override with --assembler) ----
    def illumina  = params.platform.toString().toLowerCase() == 'illumina'
    def assembler = (params.assembler ?: (illumina ? 'megahit' : 'flye')).toString().toLowerCase()
    if (illumina && assembler == 'flye')
        error "--assembler flye is long-read only; use 'megahit' or 'spades' with --platform illumina."
    if (!illumina && assembler in ['megahit', 'spades'])
        error "--assembler ${assembler} is short-read only; use 'flye' with --platform ont."

    def ch_assembled
    if (assembler == 'megahit') {
        MEGAHIT(routed.assemble)
        ch_assembled = MEGAHIT.out.assembly
    } else if (assembler == 'spades') {
        SPADES(routed.assemble)
        ch_assembled = SPADES.out.assembly
    } else {
        FLYE(routed.assemble)
        ch_assembled = FLYE.out.assembly
    }

    // short-read assemblers already produce consensus-quality contigs; racon /
    // medaka polishing only applies to the long-read path
    def ch_final
    if (illumina) {
        ch_final = ch_assembled.map { meta, bin, asm, reads -> tuple(meta, bin, asm) }
    } else {
        POLISH(ch_assembled)
        ch_final = POLISH.out.assembly
    }

    ASSEMBLY_STATS(ch_final)
    CLASSIFY_CONTIGS(ch_final.map { meta, bin, asm -> tuple(meta, bin, asm) }, ch_db)

    def all_stats = ASSEMBLY_STATS.out.stats.collectFile(
        name: 'assembly_stats.tsv', keepHeader: true, skip: 1,
        storeDir: "${params.outdir}/assembly")

    // ---- reference-consensus path ----
    if (params.reference_consensus) {
        def refcons_in = routed.refcons
            .map { meta, bin_label, taxid, reads -> tuple(taxid.toString(), meta, bin_label, reads) }
            .combine(ch_refs, by: 0)
            .map { taxid, meta, bin_label, reads, ref -> tuple(meta, bin_label, reads, ref) }
        REFERENCE_CONSENSUS(refcons_in)
    }

    emit:
    assemblies = ch_final
    stats      = all_stats
}
