version 1.0

## ont-metabuli — Terra/WDL port
##
## ONT metagenomics: QC -> host removal -> Metabuli classification ->
## significance statistics -> Sankey + HTML summary, plus an optional assembly
## module operating on Metabuli-binned reads.
##
## THIS FILE IS GENERATED.  Edit ont_metabuli.wdl.tmpl, then run:
##     python3 wdl/build_wdl.py
## The analysis scripts are embedded verbatim from bin/*.py so the Terra
## workflow needs no auxiliary file staging.
##
## Terra-specific design notes
## ---------------------------
## * Classification is BATCHED. The Metabuli database tarball (~74 GB for
##   refseq_standard) is localized ONCE and every sample is classified inside
##   that single task. Scattering classification per sample would re-localize
##   the database N times, which dominates both runtime and cost.
## * Read binning (`metabuli extract`) also needs the database, so it runs in
##   that same task. The taxa to extract must therefore be known up front:
##   the Nextflow `--assemble_taxa significant` mode depends on statistics
##   computed later and is deliberately NOT ported. List explicit taxids.
## * Containers are linux/amd64, so unlike the macOS/arm64 setup every tool
##   runs natively, and a high-memory VM lets Metabuli hold its index in RAM
##   (keep metabuli_max_ram_gb comfortably below classify_memory_gb).
## * All container images below were pulled and their tools verified.

workflow ont_metabuli {
    meta {
        description: "ONT metagenomics with Metabuli: QC, host removal, classification, significance statistics, Sankey, optional assembly."
        allowNestedInputs: true
    }

    parameter_meta {
        platform:        "'ont' (long read) or 'illumina' (paired-end). Selects the whole tool chain: chopper vs fastp, minimap2 map-ont vs sr, Metabuli --seq-mode 3 vs 2, Flye vs MEGAHIT/SPAdes. Set per RUN, not per sample."
        sample_ids:      "One entry per sample. From a Terra sample_set: this.samples.sample_id"
        fastqs:          "One FASTQ (.fastq.gz) per sample, same order as sample_ids. For platform='illumina' these are the R1 files."
        fastqs_2:        "R2 FASTQ per sample, same order as fastqs. REQUIRED when platform='illumina'; leave empty for 'ont'."
        sample_types:    "Optional per-sample type. Any of control/negative_control/ntc/blank/neg marks a negative control and switches the statistics to contaminant-aware mode."
        batches:         "Optional per-sample batch. Negative controls are matched to samples of the same batch."
        metabuli_db_tar: "Metabuli database tarball in GCS. Stage it once with stage_database.wdl."
        host_fasta:      "Host genome FASTA(.gz) to filter out, e.g. T2T-CHM13v2. Ignored when host_mmi is supplied."
        host_mmi:        "Pre-built minimap2 index of the host; skips index building. Must match the platform preset (map-ont / sr)."
        min_qual:        "ONT only. Minimum mean read quality (chopper)."
        min_len:         "ONT only. Minimum read length in bp (chopper)."
        sr_min_qual:     "Illumina only. Per-base quality threshold (fastp)."
        sr_min_len:      "Illumina only. Minimum read length after trimming (fastp)."
        adapter_fasta:   "Illumina only, optional. Adapter FASTA; fastp auto-detects for paired-end when omitted."
        assembler:       "Optional override. Defaults to flye for 'ont' and megahit for 'illumina'; 'spades' runs metaSPAdes. Long- and short-read assemblers are not interchangeable across platforms."
        flye_read_type:  "ONT only. nano-hq | nano-raw | nano-corr."
        min_contig_len:  "Illumina only. Contigs shorter than this are dropped."
        assemble_taxa:   "Taxids to assemble de novo, e.g. ['562','10239']."
        consensus_taxids: "Taxids for reference-based consensus; pair positionally with consensus_references."
    }

    input {
        # ---- sequencing platform ----
        # "ont"      : long reads  -> chopper, minimap2 map-ont, --seq-mode 3, Flye
        # "illumina" : paired-end  -> fastp,   minimap2 sr,      --seq-mode 2, MEGAHIT/SPAdes
        # Set per RUN: the statistics compare samples against each other, and
        # ONT vs Illumina differ enough in bias and depth that mixing them in one
        # significance test would be confounded. Run each platform separately.
        String         platform = "ont"

        # ---- samples (on Terra, run against a sample_set) ----
        Array[String]  sample_ids
        Array[File]    fastqs
        Array[File]?   fastqs_2      # R2 files; required when platform = "illumina"
        Array[String]? sample_types
        Array[String]? batches

        # ---- reference data ----
        File           metabuli_db_tar
        File?          host_fasta
        File?          host_mmi

        # ---- QC: long read (chopper) ----
        Boolean        skip_qc             = false
        Int            min_qual            = 10
        Int            min_len             = 500

        # ---- QC: short read (fastp) ----
        Int            sr_min_qual         = 15
        Int            sr_min_len          = 50
        File?          adapter_fasta       # optional; fastp auto-detects for PE otherwise

        # ---- host removal ----
        Boolean        skip_host_removal   = false

        # ---- classification ----
        Int            metabuli_max_ram_gb = 100
        Int            metabuli_min_score  = 0
        String         metabuli_extra      = ""

        # ---- statistics ----
        String         stats_rank          = "species"
        Int            stats_min_reads     = 10
        Float          stats_min_abundance = 0.0001
        Float          stats_alpha         = 0.05
        Float          stats_min_fold      = 5.0
        Float          stats_bg_rate       = 0.00001

        # ---- Sankey ----
        String         sankey_ranks        = "D,P,C,O,F,G,S"
        Int            sankey_top_n        = 15
        Int            sankey_min_reads    = 5

        # ---- optional assembly (all OFF by default) ----
        Boolean        assemble_unclassified = false
        Array[String]  assemble_taxa         = []
        Boolean        classify_contigs      = false
        String?        assembler             # auto: ont -> flye, illumina -> megahit; or "spades"
        String         flye_read_type        = "nano-hq"
        Int            min_bin_reads         = 50
        Int            min_contig_len        = 500    # short-read assemblers

        # ---- optional reference-based consensus (positional pairs) ----
        Array[String]  consensus_taxids      = []
        Array[File]    consensus_references  = []

        # ---- resources ----
        Int            classify_cpu        = 16
        Int            classify_memory_gb  = 120
        Int            classify_disk_gb    = 300
        String         classify_disk_type  = "SSD"
        Int            host_cpu            = 8
        Int            host_memory_gb      = 32
        Int            preemptible_tries   = 1
    }

    Int n = length(sample_ids)
    Array[File] no_files = []
    Boolean is_illumina = platform == "illumina"

    # ---- per-sample metadata (defaults when the optional arrays are absent) ----
    scatter (i in range(n)) {
        String s_type  = if defined(sample_types) then select_first([sample_types])[i] else "sample"
        String s_batch = if defined(batches)      then select_first([batches])[i]      else "all"
    }
    File metadata_tsv = write_tsv(
        flatten([[["sample", "sample_type", "batch"]], transpose([sample_ids, s_type, s_batch])])
    )

    # ---- 1. QC (platform-specific) ----
    # Reads travel as Array[File]: one entry for ONT, two (R1,R2) for Illumina.
    scatter (i in range(n)) {
        Array[File] raw_reads = if is_illumina
                                then [fastqs[i], select_first([fastqs_2])[i]]
                                else [fastqs[i]]

        if (!skip_qc && is_illumina) {
            call Fastp {
                input:
                    sample_id     = sample_ids[i],
                    reads         = raw_reads,
                    min_qual      = sr_min_qual,
                    min_len       = sr_min_len,
                    adapter_fasta = adapter_fasta,
                    preemptible_tries = preemptible_tries
            }
        }
        if (!skip_qc && !is_illumina) {
            call Chopper {
                input:
                    sample_id = sample_ids[i],
                    reads     = fastqs[i],
                    min_qual  = min_qual,
                    min_len   = min_len,
                    preemptible_tries = preemptible_tries
            }
        }
        Array[File] qcd_reads = select_first([Fastp.filtered, Chopper.filtered, raw_reads])
        File?       qc_json   = if is_illumina then Fastp.report else Chopper.report
    }

    # ---- 2. host removal (index built once, filtering scattered) ----
    if (!skip_host_removal) {
        if (!defined(host_mmi)) {
            call Minimap2Index {
                input:
                    host_fasta = select_first([host_fasta]),
                    preset     = if is_illumina then "sr" else "map-ont",
                    preemptible_tries = preemptible_tries
            }
        }
        File mmi = select_first([host_mmi, Minimap2Index.index])

        scatter (i in range(n)) {
            call HostRemoval {
                input:
                    sample_id = sample_ids[i],
                    reads     = qcd_reads[i],
                    paired    = is_illumina,
                    index     = mmi,
                    cpu       = host_cpu,
                    memory_gb = host_memory_gb,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    Array[Array[File]] clean_reads = select_first([HostRemoval.filtered, qcd_reads])
    Array[File] host_flagstat_files = select_first([HostRemoval.flagstat, no_files])

    # split into parallel R1 / R2 arrays for the batched classifier
    scatter (cr in clean_reads) { File clean_r1 = cr[0] }
    if (is_illumina) { scatter (cr in clean_reads) { File clean_r2_i = cr[1] } }
    Array[File] clean_r2 = select_first([clean_r2_i, no_files])

    # ---- 3. classification + read binning (one database localization) ----
    call MetabuliClassifyAndExtract {
        input:
            sample_ids           = sample_ids,
            reads                = clean_r1,
            reads_2              = clean_r2,
            paired               = is_illumina,
            db_tar               = metabuli_db_tar,
            max_ram_gb           = metabuli_max_ram_gb,
            min_score            = metabuli_min_score,
            extra                = metabuli_extra,
            extract_unclassified = assemble_unclassified,
            extract_taxa_denovo  = assemble_taxa,
            extract_taxa_ref     = consensus_taxids,
            cpu                  = classify_cpu,
            memory_gb            = classify_memory_gb,
            disk_gb              = classify_disk_gb,
            disk_type            = classify_disk_type
    }

    # ---- 4. combined abundance tables ----
    call CombineReports {
        input:
            reports  = MetabuliClassifyAndExtract.reports,
            metadata = metadata_tsv,
            rank     = stats_rank,
            preemptible_tries = preemptible_tries
    }

    # ---- 5. significance / decontamination statistics ----
    call DecontamStats {
        input:
            abundance     = CombineReports.long_table,
            metadata      = metadata_tsv,
            rank          = stats_rank,
            min_reads     = stats_min_reads,
            min_abundance = stats_min_abundance,
            alpha         = stats_alpha,
            min_fold      = stats_min_fold,
            bg_rate       = stats_bg_rate,
            preemptible_tries = preemptible_tries
    }

    # ---- 6. per-sample Sankey ----
    scatter (report in MetabuliClassifyAndExtract.reports) {
        call Sankey {
            input:
                report    = report,
                sample_id = basename(report, "_report.tsv"),
                flags     = DecontamStats.flags,
                ranks     = sankey_ranks,
                top_n     = sankey_top_n,
                min_reads = sankey_min_reads,
                preemptible_tries = preemptible_tries
        }
    }

    # ---- 7. optional de novo assembly (bins not labelled ref_*) ----
    # Bins are emitted as <sample>__<label>_R1.fastq (plus _R2 when paired), so
    # the R1 and R2 globs sort into the same stem order and pair by index.
    String asm = select_first([assembler, if is_illumina then "megahit" else "flye"])

    scatter (bi in range(length(MetabuliClassifyAndExtract.bins_r1))) {
        File   dn_r1    = MetabuliClassifyAndExtract.bins_r1[bi]
        String dn_stem  = basename(dn_r1, "_R1.fastq")
        String dn_label = sub(dn_stem, "^.*__", "")
        String dn_sample = sub(dn_stem, "__.*$", "")
        Array[File] dn_reads = if is_illumina
                               then [dn_r1, MetabuliClassifyAndExtract.bins_r2[bi]]
                               else [dn_r1]

        # a de novo bin is any bin whose label does not start with ref_
        if (sub(dn_label, "^ref_", "") == dn_label) {
            if (asm == "flye") {
                call Flye {
                    input:
                        sample_id     = dn_sample,
                        bin_label     = dn_label,
                        reads         = dn_r1,
                        read_type     = flye_read_type,
                        min_bin_reads = min_bin_reads,
                        preemptible_tries = preemptible_tries
                }
            }
            if (asm != "flye") {
                call ShortReadAssembly {
                    input:
                        sample_id      = dn_sample,
                        bin_label      = dn_label,
                        reads          = dn_reads,
                        assembler      = asm,
                        min_bin_reads  = min_bin_reads,
                        min_contig_len = min_contig_len,
                        preemptible_tries = preemptible_tries
                }
            }
            File bin_assembly = select_first([Flye.assembly, ShortReadAssembly.assembly])
            call AssemblyStats {
                input:
                    sample_id = dn_sample,
                    bin_label = dn_label,
                    assembly  = bin_assembly,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    Array[File] assemblies = select_all(bin_assembly)

    if (classify_contigs) {
        call ClassifyContigs {
            input:
                assemblies = assemblies,
                db_tar     = metabuli_db_tar,
                max_ram_gb = metabuli_max_ram_gb,
                cpu        = classify_cpu,
                memory_gb  = classify_memory_gb,
                disk_gb    = classify_disk_gb,
                disk_type  = classify_disk_type
        }
    }

    # ---- 8. optional reference-based consensus ----
    # Single-level scatter over the ref_* bins. The reference for a bin is
    # chosen inside the task from the (taxid, reference) pairs, which avoids a
    # nested scatter — Cromwell cannot flatten optionals out of one.
    scatter (ri in range(length(MetabuliClassifyAndExtract.bins_r1))) {
        File   r_r1    = MetabuliClassifyAndExtract.bins_r1[ri]
        String r_stem  = basename(r_r1, "_R1.fastq")
        String r_label = sub(r_stem, "^.*__", "")
        Array[File] r_reads = if is_illumina
                              then [r_r1, MetabuliClassifyAndExtract.bins_r2[ri]]
                              else [r_r1]
        if (sub(r_label, "^ref_", "") != r_label) {
            call MapToReference {
                input:
                    sample_id  = sub(r_stem, "__.*$", ""),
                    bin_label  = r_label,
                    taxid      = sub(r_label, "^ref_", ""),
                    reads      = r_reads,
                    paired     = is_illumina,
                    ref_taxids = consensus_taxids,
                    references = consensus_references,
                    cpu        = host_cpu,
                    memory_gb  = host_memory_gb,
                    preemptible_tries = preemptible_tries
            }
            call CallVariants {
                input:
                    sample_id     = sub(r_stem, "__.*$", ""),
                    bin_label     = r_label,
                    bam           = MapToReference.bam,
                    bam_index     = MapToReference.bam_index,
                    reference     = MapToReference.prepared_reference,
                    reference_fai = MapToReference.prepared_reference_fai,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    # ---- 9. HTML summary ----
    call Summary {
        input:
            sankeys      = Sankey.html,
            kronas       = MetabuliClassifyAndExtract.kronas,
            significant  = DecontamStats.significant,
            summary_json = DecontamStats.summary_json,
            qc_jsons     = select_all(qc_json),
            host_stats   = host_flagstat_files,
            preemptible_tries = preemptible_tries
    }

    output {
        # classification
        Array[File] metabuli_reports         = MetabuliClassifyAndExtract.reports
        Array[File] metabuli_classifications = MetabuliClassifyAndExtract.classifications
        Array[File] krona_html               = MetabuliClassifyAndExtract.kronas

        # abundance + statistics
        File        abundance_long           = CombineReports.long_table
        File        abundance_matrix         = CombineReports.matrix
        File        read_accounting          = CombineReports.accounting
        File        taxon_stats              = DecontamStats.table
        File        significant_taxa         = DecontamStats.significant
        File        stats_summary_json       = DecontamStats.summary_json

        # visualisation
        Array[File] sankey_html              = Sankey.html
        File        summary_html             = Summary.html
        File        run_summary              = Summary.table

        # QC / host removal
        Array[File] qc_reports               = select_all(qc_json)
        Array[File] host_flagstats           = host_flagstat_files

        # optional assembly
        Array[File] assembly_fasta           = assemblies
        Array[File] assembly_stats           = select_all(AssemblyStats.stats)
        File?       contig_classification    = ClassifyContigs.merged_report

        # optional reference-based consensus
        Array[File] consensus_fasta          = select_all(MapToReference.consensus)
        Array[File] consensus_coverage       = select_all(MapToReference.coverage)
        Array[File] consensus_vcf            = select_all(CallVariants.vcf)
    }
}

# ============================================================================
#  Tasks
# ============================================================================

task Chopper {
    input {
        String sample_id
        File   reads
        Int    min_qual
        Int    min_len
        Int    preemptible_tries
        Int    cpu = 4
        Int    memory_gb = 8
    }
    Int disk_gb = ceil(size(reads, "GB") * 4) + 20

    command <<<
        set -euo pipefail
        chopper --quality ~{min_qual} --minlength ~{min_len} --threads ~{cpu} \
            --input ~{reads} 2> chopper.log | gzip > ~{sample_id}.filtered.fastq.gz

        # distil chopper's "Kept X reads out of Y reads" line into a small JSON
        kept=$(awk '/Kept/{print $2}'  chopper.log); kept=${kept:-0}
        total=$(awk '/Kept/{print $6}' chopper.log); total=${total:-0}
        printf '{"sample":"%s","input_reads":%s,"reads":%s,"min_qual":%s,"min_len":%s}\n' \
            "~{sample_id}" "$total" "$kept" "~{min_qual}" "~{min_len}" > ~{sample_id}.qc.json
    >>>

    output {
        Array[File] filtered = ["~{sample_id}.filtered.fastq.gz"]
        File        report   = "~{sample_id}.qc.json"
    }
    runtime {
        docker: "quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible_tries
    }
}

task Fastp {
    input {
        String      sample_id
        Array[File] reads          # [R1, R2]
        Int         min_qual
        Int         min_len
        File?       adapter_fasta
        Int         preemptible_tries
        Int         cpu = 4
        Int         memory_gb = 8
    }
    Int disk_gb = ceil(size(reads, "GB") * 4) + 20

    command <<<
        set -euo pipefail
        # Paired-end adapter trimming + quality filtering. fastp auto-detects
        # adapters for PE data unless an adapter FASTA is given.
        fastp \
            --in1 ~{reads[0]} --in2 ~{reads[1]} \
            --out1 ~{sample_id}_R1.filtered.fastq.gz \
            --out2 ~{sample_id}_R2.filtered.fastq.gz \
            ~{if defined(adapter_fasta) then "--adapter_fasta " + adapter_fasta else "--detect_adapter_for_pe"} \
            --qualified_quality_phred ~{min_qual} \
            --length_required ~{min_len} \
            --thread ~{cpu} \
            --json ~{sample_id}.fastp.json \
            --html ~{sample_id}.fastp.html

        # Normalise fastp's JSON into the small schema the summary step reads.
        # awk keeps this independent of what interpreters the image ships.
        get_total () {
            awk -v sec="\"$1\"" 'index($0, sec){f=1} f && /"total_reads"/ {gsub(/[^0-9]/,""); print; exit}' \
                ~{sample_id}.fastp.json
        }
        before=$(get_total before_filtering); before=${before:-0}
        after=$(get_total after_filtering);   after=${after:-0}
        printf '{"sample":"%s","input_reads":%s,"reads":%s,"min_qual":%s,"min_len":%s}\n' \
            "~{sample_id}" "$before" "$after" "~{min_qual}" "~{min_len}" > ~{sample_id}.qc.json
    >>>

    output {
        Array[File] filtered = ["~{sample_id}_R1.filtered.fastq.gz", "~{sample_id}_R2.filtered.fastq.gz"]
        File        report   = "~{sample_id}.qc.json"
        File        html     = "~{sample_id}.fastp.html"
    }
    runtime {
        docker: "quay.io/biocontainers/fastp:1.3.6--h43da1c4_0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible_tries
    }
}

task Minimap2Index {
    input {
        File   host_fasta
        String preset          # map-ont (long read) | sr (short read)
        Int    preemptible_tries
        Int    cpu = 8
        Int    memory_gb = 32
    }
    Int disk_gb = ceil(size(host_fasta, "GB") * 6) + 30

    command <<<
        set -euo pipefail
        # the index must be built with the same preset used for mapping
        minimap2 -x ~{preset} -d host.mmi -t ~{cpu} ~{host_fasta}
    >>>

    output { File index = "host.mmi" }
    runtime {
        docker: "quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible_tries
    }
}

task HostRemoval {
    input {
        String      sample_id
        Array[File] reads         # [R1] or [R1, R2]
        Boolean     paired
        File        index
        Int         cpu
        Int         memory_gb
        Int         preemptible_tries
    }
    Int disk_gb = ceil(size(reads, "GB") * 8 + size(index, "GB") * 2) + 40

    command <<<
        set -euo pipefail
        # Keep only reads that do NOT map to the host. Streamed through a
        # temporary BAM so a full host alignment is never materialised.
        minimap2 -ax ~{if paired then "sr" else "map-ont"} -t ~{cpu} \
            ~{index} ~{sep=' ' reads} 2> minimap2.log \
            | samtools view -@ ~{cpu} -b -o tmp.bam -

        samtools flagstat -@ ~{cpu} tmp.bam > ~{sample_id}.host.flagstat

        if ~{if paired then "true" else "false"}; then
            # keep pairs where BOTH mates are unmapped (-f 12); collate restores
            # mate adjacency so samtools fastq can emit a proper R1/R2 set
            samtools view -@ ~{cpu} -b -f 12 -F 256 tmp.bam \
                | samtools collate -@ ~{cpu} -u -O - \
                | samtools fastq -@ ~{cpu} -n \
                    -1 ~{sample_id}_R1.hostremoved.fastq.gz \
                    -2 ~{sample_id}_R2.hostremoved.fastq.gz \
                    -0 /dev/null -s /dev/null -
        else
            samtools fastq -@ ~{cpu} -n -f 4 tmp.bam | gzip > ~{sample_id}_R1.hostremoved.fastq.gz
        fi
        rm -f tmp.bam
    >>>

    output {
        Array[File] filtered = if paired
            then ["~{sample_id}_R1.hostremoved.fastq.gz", "~{sample_id}_R2.hostremoved.fastq.gz"]
            else ["~{sample_id}_R1.hostremoved.fastq.gz"]
        File flagstat = "~{sample_id}.host.flagstat"
    }
    runtime {
        docker: "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible_tries
    }
}

task MetabuliClassifyAndExtract {
    input {
        Array[String] sample_ids
        Array[File]   reads          # R1 (or the single ONT file) per sample
        Array[File]   reads_2        # R2 per sample; empty for ONT
        Boolean       paired
        File          db_tar
        Int           max_ram_gb
        Int           min_score
        String        extra
        Boolean       extract_unclassified
        Array[String] extract_taxa_denovo
        Array[String] extract_taxa_ref
        Int           cpu
        Int           memory_gb
        Int           disk_gb
        String        disk_type
    }

    command <<<
        set -euo pipefail

        echo "== unpacking Metabuli database =="
        mkdir -p dbroot out reports bins
        tar -xzf ~{db_tar} -C dbroot
        # the tarball may nest the database one or more levels deep, so locate
        # it by its parameter file instead of assuming a layout
        DBDIR=$(dirname "$(find dbroot -maxdepth 4 -name db.parameters | head -1)")
        if [ -z "$DBDIR" ] || [ ! -d "$DBDIR" ]; then
            echo "ERROR: no Metabuli database found inside the tarball" >&2
            find dbroot -maxdepth 3 >&2
            exit 1
        fi
        echo "database: $DBDIR"
        cat "$DBDIR/db.parameters" || true

        ids=(~{sep=' ' sample_ids})
        rds=(~{sep=' ' reads})
        rd2=(~{sep=' ' reads_2})
        PAIRED=~{if paired then "1" else "0"}
        SEQMODE=~{if paired then "2" else "3"}

        # Bin one clade (or the unclassified pile) out of a sample's reads.
        # metabuli names its output after the INPUT basenames, so inputs are
        # staged under canonical names to make the results predictable.
        extract_bin () {   # $1 sample, $2 classifications, $3 taxid, $4 label
            rm -rf ext && mkdir -p ext
            if [ "$PAIRED" = "1" ]; then
                metabuli extract q_R1.fastq.gz q_R2.fastq.gz "$2" "$DBDIR" \
                    --tax-id "$3" --seq-mode 2 --extract-format 2 --outdir ext || true
                if [ -s "ext/q_R1_$3.fq" ] && [ -s "ext/q_R2_$3.fq" ]; then
                    mv "ext/q_R1_$3.fq" "bins/$1__$4_R1.fastq"
                    mv "ext/q_R2_$3.fq" "bins/$1__$4_R2.fastq"
                    echo "   bin $4: $(( $(wc -l < "bins/$1__$4_R1.fastq") / 4 )) pairs"
                else
                    echo "   bin $4: empty"
                fi
            else
                metabuli extract q_R1.fastq.gz "$2" "$DBDIR" \
                    --tax-id "$3" --seq-mode 3 --extract-format 2 --outdir ext || true
                found=$(ls ext/*.fq 2>/dev/null | head -1 || true)
                if [ -n "$found" ]; then
                    mv "$found" "bins/$1__$4_R1.fastq"
                    echo "   bin $4: $(( $(wc -l < "bins/$1__$4_R1.fastq") / 4 )) reads"
                else
                    echo "   bin $4: empty"
                fi
            fi
        }

        for i in "${!ids[@]}"; do
            sid="${ids[$i]}"
            echo "== classifying $sid =="

            # canonical staging (also what extract_bin reads)
            rm -f q_R1.fastq.gz q_R2.fastq.gz
            ln -s "${rds[$i]}" q_R1.fastq.gz
            QUERY="q_R1.fastq.gz"
            if [ "$PAIRED" = "1" ]; then
                ln -s "${rd2[$i]}" q_R2.fastq.gz
                QUERY="q_R1.fastq.gz q_R2.fastq.gz"
            fi

            metabuli classify \
                --seq-mode "$SEQMODE" \
                ~{if paired then "--precise 1" else ""} \
                --threads ~{cpu} \
                --max-ram ~{max_ram_gb} \
                ~{if min_score > 0 then "--min-score " + min_score else ""} \
                ~{extra} \
                $QUERY "$DBDIR" out "$sid"

            cp "out/${sid}_report.tsv" "reports/${sid}_report.tsv"
            if [ -f "out/${sid}_classifications.tsv" ]; then
                cp "out/${sid}_classifications.tsv" "reports/${sid}_classifications.tsv"
            fi
            if [ -f "out/${sid}_krona.html" ]; then
                cp "out/${sid}_krona.html" "reports/${sid}_krona.html"
            fi

            # ---- read binning for the optional assembly / consensus modules ----
            # performed here because `metabuli extract` needs the database, which
            # is already unpacked in this task
            cls="out/${sid}_classifications.tsv"
            if [ -f "$cls" ]; then
                if ~{if extract_unclassified then "true" else "false"}; then
                    extract_bin "$sid" "$cls" -1 unclassified
                fi
                for t in ~{sep=' ' extract_taxa_denovo}; do
                    extract_bin "$sid" "$cls" "$t" "taxon_$t"
                done
                for t in ~{sep=' ' extract_taxa_ref}; do
                    extract_bin "$sid" "$cls" "$t" "ref_$t"
                done
            fi
        done
    >>>

    output {
        Array[File] reports         = glob("reports/*_report.tsv")
        Array[File] classifications = glob("reports/*_classifications.tsv")
        Array[File] kronas          = glob("reports/*_krona.html")
        # R1 and R2 globs share the same stems, so they sort into matching order
        Array[File] bins_r1         = glob("bins/*_R1.fastq")
        Array[File] bins_r2         = glob("bins/*_R2.fastq")
    }
    runtime {
        docker: "quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " " + disk_type
        preemptible: 0
        maxRetries: 1
    }
}

task CombineReports {
    input {
        Array[File] reports
        File        metadata
        String      rank
        Int         preemptible_tries
    }

    command <<<
        set -euo pipefail
        pip install --no-cache-dir --quiet 'pandas==2.2.3' 'numpy==2.1.3' 1>&2

cat > combine_reports.py <<'PYEOF'
#!/usr/bin/env python3
"""
Combine per-sample Metabuli/Kraken2-style report.tsv files into tidy abundance
tables. Produces a long table (one row per sample x taxon at the target rank),
a wide count matrix, and a per-sample read-accounting table.

Report columns (tab-separated, no header):
    pct  clade_reads  taxon_reads  rank_code  taxid  name(indented 2sp/level)
"""
import argparse
import os
import sys

RANK_LETTER = {'species': 'S', 'genus': 'G', 'family': 'F', 'order': 'O',
               'class': 'C', 'phylum': 'P', 'domain': 'D', 'kingdom': 'K'}

# Metabuli v1.2.0 reports use full rank names ("species", "domain", ...);
# older/Kraken2 style uses single letters ("S", "D", ...). Normalise both.
_FULL2LETTER = {'superkingdom': 'D', 'domain': 'D', 'kingdom': 'K', 'phylum': 'P',
                'class': 'C', 'order': 'O', 'family': 'F', 'genus': 'G', 'species': 'S'}


def norm_rank(rank):
    """Return canonical single-letter code, or None for sub-ranks/no-rank."""
    r = (rank or '').strip()
    if len(r) == 1 and r.upper() in 'DKPCOFGS':
        return r.upper()
    return _FULL2LETTER.get(r.lower())


def parse_report(path):
    """Return (rows, totals) for one report.
    rows: list of dicts with taxid, rank, name, taxon_reads, clade_reads, depth.
    totals: dict with classified / unclassified / total.
    """
    rows = []
    classified = unclassified = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            try:
                clade_reads = int(parts[1])
                taxon_reads = int(parts[2])
            except ValueError:
                continue  # header or malformed
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name_raw = parts[5]
            depth = (len(name_raw) - len(name_raw.lstrip(' '))) // 2
            name = name_raw.strip()
            if rank == 'U' or taxid == '0':
                unclassified += clade_reads
                continue
            if rank == 'R' or taxid == '1':
                classified = max(classified, clade_reads)
            rows.append({'taxid': taxid, 'rank': rank, 'name': name,
                         'taxon_reads': taxon_reads, 'clade_reads': clade_reads,
                         'depth': depth})
    if classified == 0:  # no explicit root row: sum top-level clades
        classified = sum(r['clade_reads'] for r in rows if r['depth'] <= 1)
    totals = {'classified': classified, 'unclassified': unclassified,
              'total': classified + unclassified}
    return rows, totals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reports', nargs='+', required=True)
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--rank', default='species')
    ap.add_argument('--out-prefix', default='combined')
    args = ap.parse_args()

    target = RANK_LETTER.get(args.rank.lower(), 'S')

    long_rows = []
    accounting = []
    matrix = {}      # (taxid, name) -> {sample: clade_reads}
    samples = []

    for rep in args.reports:
        base = os.path.basename(rep)
        sample = base[:-len('_report.tsv')] if base.endswith('_report.tsv') else base.split('.')[0]
        samples.append(sample)
        rows, totals = parse_report(rep)
        tot = totals['classified'] or 1
        accounting.append((sample, totals['total'], totals['classified'],
                           totals['unclassified'],
                           round(100.0 * totals['classified'] / (totals['total'] or 1), 2)))
        for r in rows:
            # match on normalised rank so both "S" and "species" work; sub-ranks
            # (strain/serogroup/no rank) normalise to None and are excluded
            if norm_rank(r['rank']) == target:
                rel = r['clade_reads'] / tot
                long_rows.append((sample, r['taxid'], target, r['name'],
                                  r['taxon_reads'], r['clade_reads'],
                                  totals['classified'], f"{rel:.6g}"))
                key = (r['taxid'], r['name'])
                matrix.setdefault(key, {})[sample] = r['clade_reads']

    # ---- long table ----
    with open('combined_abundance_long.tsv', 'w') as fh:
        fh.write('sample\ttaxid\trank\tname\ttaxon_reads\tclade_reads\ttotal_classified\trel_abundance\n')
        for row in long_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    # ---- wide matrix ----
    samples_sorted = sorted(set(samples))
    with open(f'abundance_matrix_{args.rank.lower()}.tsv', 'w') as fh:
        fh.write('taxid\tname\t' + '\t'.join(samples_sorted) + '\n')
        for (taxid, name), d in sorted(matrix.items(), key=lambda kv: -sum(kv[1].values())):
            counts = [str(d.get(s, 0)) for s in samples_sorted]
            fh.write(f'{taxid}\t{name}\t' + '\t'.join(counts) + '\n')

    # ---- read accounting ----
    with open('read_accounting.tsv', 'w') as fh:
        fh.write('sample\ttotal_reads\tclassified\tunclassified\tpct_classified\n')
        for row in sorted(accounting):
            fh.write('\t'.join(str(x) for x in row) + '\n')

    print(f"[combine_reports] {len(samples_sorted)} samples, "
          f"{len(matrix)} {args.rank}-level taxa", file=sys.stderr)


if __name__ == '__main__':
    main()
PYEOF

        python3 combine_reports.py \
            --reports ~{sep=' ' reports} \
            --metadata ~{metadata} \
            --rank ~{rank} \
            --out-prefix combined
    >>>

    output {
        File long_table = "combined_abundance_long.tsv"
        File matrix     = glob("abundance_matrix_*.tsv")[0]
        File accounting = "read_accounting.tsv"
    }
    runtime {
        docker: "python:3.11-slim"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk 50 HDD"
        preemptible: preemptible_tries
    }
}

task DecontamStats {
    input {
        File   abundance
        File   metadata
        String rank
        Int    min_reads
        Float  min_abundance
        Float  alpha
        Float  min_fold
        Float  bg_rate
        Int    preemptible_tries
    }

    command <<<
        set -euo pipefail
        pip install --no-cache-dir --quiet 'pandas==2.2.3' 'numpy==2.1.3' 'scipy==1.14.1' 1>&2

cat > decontam_stats.py <<'PYEOF'
#!/usr/bin/env python3
"""
Significance / decontamination statistics for ONT metagenomics.

Two modes, auto-selected from the samplesheet:

  (A) NEGATIVE-CONTROL mode  — when >=1 sample is flagged as a control.
      For each biological sample and each taxon, compare the taxon's proportion
      of classified reads in the sample vs. in the (batch-matched, else pooled)
      negative control with a one-sided Fisher exact test (enriched in sample).
      A taxon is called:
        contaminant       if its proportion in the control >= in the sample
        significant       if BH q < alpha AND fold-enrichment >= min_fold
        below_threshold   if reads < min_reads or rel-abundance < min_abundance
        not_significant   otherwise
      This is the frequency-based logic behind decontam-style contaminant
      identification, adapted for a single control per batch.

  (B) NO-CONTROL mode  — when no control is present (the common ONT case).
      Significance testing without a blank is inherently limited: there is no
      empirical contamination background to test against. We therefore report a
      transparent HEURISTIC: for each taxon, a one-sided Poisson tail test of the
      observed read count against an assumed low background rate (--bg-rate x
      classified reads), plus a Wilson 95% CI on relative abundance. A taxon is:
        detected          if reads>=min_reads, rel-abundance>=min_abundance, BH q<alpha
        below_threshold   if it fails the count/abundance floor
        not_significant   otherwise
      Treat these as detection confidence flags, NOT contamination-corrected
      calls. Include a negative control whenever possible for rigorous results.

Outputs: taxon_stats.tsv, significant_taxa.tsv, sankey_flags.tsv, stats_summary.json
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, poisson

CONTROL_CALL = 'contaminant'


def bh_qvalues(pvals):
    """Benjamini-Hochberg FDR. Returns array aligned to input order."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def wilson_ci(count, total, z=1.96):
    if total == 0:
        return (0.0, 0.0)
    phat = count / total
    denom = 1 + z**2 / total
    centre = (phat + z**2 / (2 * total)) / denom
    half = (z * math.sqrt(phat * (1 - phat) / total + z**2 / (4 * total**2))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--abundance', required=True)   # combined_abundance_long.tsv
    ap.add_argument('--metadata', required=True)    # sample sample_type batch
    ap.add_argument('--rank', default='species')
    ap.add_argument('--min-reads', type=int, default=10)
    ap.add_argument('--min-abundance', type=float, default=1e-4)
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--min-fold', type=float, default=5.0)
    ap.add_argument('--bg-rate', type=float, default=1e-5)
    ap.add_argument('--out-prefix', default='.')
    args = ap.parse_args()
    os.makedirs(args.out_prefix, exist_ok=True)
    eps = 1e-9

    abund = pd.read_csv(args.abundance, sep='\t', dtype={'taxid': str})
    meta = pd.read_csv(args.metadata, sep='\t', dtype=str).fillna('all')
    meta['sample_type'] = meta['sample_type'].str.lower()
    controls = meta.loc[meta['sample_type'] == 'control', 'sample'].tolist()
    bio_samples = meta.loc[meta['sample_type'] != 'control', 'sample'].tolist()
    batch_of = dict(zip(meta['sample'], meta['batch']))

    # per-sample classified total (constant within a sample)
    totals = abund.groupby('sample')['total_classified'].max().to_dict()
    # reads[sample][taxid] = clade_reads ; names[taxid] = name
    reads = {}
    names = {}
    for _, r in abund.iterrows():
        reads.setdefault(r['sample'], {})[r['taxid']] = int(r['clade_reads'])
        names[r['taxid']] = r['name']

    mode = 'negative_control' if controls else 'no_control'
    out_rows = []
    per_sample_sig = {}

    def control_pool(sample):
        """Return (reads_dict, total) for controls matched to this sample."""
        b = batch_of.get(sample, 'all')
        matched = [c for c in controls if batch_of.get(c, 'all') == b]
        if not matched:
            matched = controls  # fall back to all controls pooled
        pooled = {}
        tot = 0
        for c in matched:
            tot += totals.get(c, 0)
            for tx, ct in reads.get(c, {}).items():
                pooled[tx] = pooled.get(tx, 0) + ct
        return pooled, tot

    for sample in bio_samples:
        A = totals.get(sample, 0)
        if A == 0:
            continue
        s_reads = reads.get(sample, {})
        taxa = list(s_reads.keys())

        if mode == 'negative_control':
            c_reads, C = control_pool(sample)
            pvals, recs = [], []
            for tx in taxa:
                a = s_reads[tx]
                c = c_reads.get(tx, 0)
                prop_s = a / A
                prop_c = (c / C) if C else 0.0
                table = [[a, max(A - a, 0)], [c, max(C - c, 0)]]
                try:
                    _, p = fisher_exact(table, alternative='greater')
                except ValueError:
                    p = 1.0
                fold = (prop_s + eps) / (prop_c + eps)
                pvals.append(p)
                recs.append((tx, a, prop_s, c, prop_c, fold, p))
            qvals = bh_qvalues(pvals)
            nsig = 0
            for (tx, a, prop_s, c, prop_c, fold, p), q in zip(recs, qvals):
                if a < args.min_reads or prop_s < args.min_abundance:
                    call = 'below_threshold'
                elif prop_c >= prop_s:
                    call = CONTROL_CALL
                elif q < args.alpha and fold >= args.min_fold:
                    call = 'significant'
                else:
                    call = 'not_significant'
                if call == 'significant':
                    nsig += 1
                lo, hi = wilson_ci(a, A)
                out_rows.append([sample, tx, names.get(tx, tx), args.rank, a, A,
                                 f'{prop_s:.6g}', f'{lo:.4g}', f'{hi:.4g}',
                                 c, f'{prop_c:.6g}', f'{fold:.4g}',
                                 f'{p:.4g}', f'{q:.4g}', call, mode])
            per_sample_sig[sample] = nsig

        else:  # no_control
            pvals, recs = [], []
            for tx in taxa:
                a = s_reads[tx]
                prop_s = a / A
                lam = max(args.bg_rate * A, eps)
                p = float(poisson.sf(a - 1, lam))  # P(X >= a)
                pvals.append(p)
                recs.append((tx, a, prop_s, p))
            qvals = bh_qvalues(pvals)
            nsig = 0
            for (tx, a, prop_s, p), q in zip(recs, qvals):
                if a < args.min_reads or prop_s < args.min_abundance:
                    call = 'below_threshold'
                elif q < args.alpha:
                    call = 'detected'
                else:
                    call = 'not_significant'
                if call == 'detected':
                    nsig += 1
                lo, hi = wilson_ci(a, A)
                out_rows.append([sample, tx, names.get(tx, tx), args.rank, a, A,
                                 f'{prop_s:.6g}', f'{lo:.4g}', f'{hi:.4g}',
                                 '', '', '', f'{p:.4g}', f'{q:.4g}', call, mode])
            per_sample_sig[sample] = nsig

    cols = ['sample', 'taxid', 'name', 'rank', 'reads', 'total_classified',
            'rel_abundance', 'ci_low', 'ci_high', 'control_reads',
            'control_rel_abundance', 'fold_enrichment', 'p_value', 'q_value',
            'call', 'mode']
    df = pd.DataFrame(out_rows, columns=cols)
    df = df.sort_values(['sample', 'reads'], ascending=[True, False])

    p = lambda f: os.path.join(args.out_prefix, f)
    df.to_csv(p('taxon_stats.tsv'), sep='\t', index=False)

    sig_calls = {'significant', 'detected'}
    df[df['call'].isin(sig_calls)].to_csv(p('significant_taxa.tsv'), sep='\t', index=False)

    df[['sample', 'taxid', 'call']].to_csv(p('sankey_flags.tsv'), sep='\t', index=False)

    summary = {
        'mode': mode,
        'n_biological_samples': len(bio_samples),
        'n_controls': len(controls),
        'controls': controls,
        'rank': args.rank,
        'params': {'min_reads': args.min_reads, 'min_abundance': args.min_abundance,
                   'alpha': args.alpha, 'min_fold': args.min_fold,
                   'bg_rate': args.bg_rate},
        'significant_per_sample': per_sample_sig,
        'n_significant_total': int(sum(per_sample_sig.values())),
        'caveat': ('No negative control supplied: calls are Poisson-based '
                   'detection-confidence heuristics, not contamination-corrected.'
                   if mode == 'no_control' else
                   'Contaminant calls are frequency-based vs. batch-matched controls.')
    }
    with open(p('stats_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2)

    print(f"[decontam_stats] mode={mode} controls={len(controls)} "
          f"sig_total={summary['n_significant_total']}", file=sys.stderr)


if __name__ == '__main__':
    main()
PYEOF

        python3 decontam_stats.py \
            --abundance ~{abundance} \
            --metadata ~{metadata} \
            --rank ~{rank} \
            --min-reads ~{min_reads} \
            --min-abundance ~{min_abundance} \
            --alpha ~{alpha} \
            --min-fold ~{min_fold} \
            --bg-rate ~{bg_rate} \
            --out-prefix .
    >>>

    output {
        File table        = "taxon_stats.tsv"
        File significant  = "significant_taxa.tsv"
        File flags        = "sankey_flags.tsv"
        File summary_json = "stats_summary.json"
    }
    runtime {
        docker: "python:3.11-slim"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk 50 HDD"
        preemptible: preemptible_tries
    }
}

task Sankey {
    input {
        File   report
        String sample_id
        File   flags
        String ranks
        Int    top_n
        Int    min_reads
        Int    preemptible_tries
    }

    command <<<
        set -euo pipefail
        pip install --no-cache-dir --quiet 'plotly==5.24.1' 1>&2

cat > make_sankey.py <<'PYEOF'
#!/usr/bin/env python3
"""
Build a self-contained interactive taxonomic Sankey (Domain -> ... -> Species)
from one Metabuli/Kraken2-style report.tsv, using plotly.

Flow width = clade reads. Nodes are laid out in rank columns. Species nodes are
recoloured by their decontam call (from sankey_flags.tsv) when available:
significant/detected = green, contaminant = grey, below/ns = muted.

The output HTML embeds plotly.js (include_plotlyjs=True) so it works offline.
"""
import argparse
import os
import sys

# canonical single-letter rank codes we may place as columns
RANK_NAMES = {'D': 'Domain', 'K': 'Kingdom', 'P': 'Phylum', 'C': 'Class',
              'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'}

# Metabuli v1.2.0 uses full rank names; older format uses single letters.
_FULL2LETTER = {'superkingdom': 'D', 'domain': 'D', 'kingdom': 'K', 'phylum': 'P',
                'class': 'C', 'order': 'O', 'family': 'F', 'genus': 'G', 'species': 'S'}


def norm_rank(rank):
    r = (rank or '').strip()
    if len(r) == 1 and r.upper() in 'DKPCOFGS':
        return r.upper()
    return _FULL2LETTER.get(r.lower())

DOMAIN_COLORS = {
    'Bacteria': '#4C78A8', 'Viruses': '#E45756', 'Archaea': '#B279A2',
    'Eukaryota': '#F58518', 'Fungi': '#72B7B2', 'unknown': '#9D9D9D',
}
CALL_COLORS = {
    'significant': '#2CA02C', 'detected': '#2CA02C',
    'contaminant': '#8C8C8C', 'below_threshold': '#D9D9D9',
    'not_significant': '#BFBFBF',
}


def parse_tree(path):
    """Parse report into ordered nodes with parent pointers (via depth stack)."""
    nodes = []
    stack = {}  # depth -> node index
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            try:
                clade_reads = int(parts[1])
            except ValueError:
                continue
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name_raw = parts[5]
            depth = (len(name_raw) - len(name_raw.lstrip(' '))) // 2
            name = name_raw.strip()
            if taxid == '0':  # unclassified
                continue
            parent = None
            for d in range(depth - 1, -1, -1):
                if d in stack:
                    parent = stack[d]
                    break
            idx = len(nodes)
            nodes.append({'taxid': taxid, 'rank': rank, 'name': name,
                          'clade_reads': clade_reads, 'depth': depth,
                          'parent': parent})
            stack[depth] = idx
            # drop deeper stack entries
            for d in list(stack):
                if d > depth:
                    del stack[d]
    return nodes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--report', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--flags', default=None)
    ap.add_argument('--ranks', default='D,P,C,O,F,G,S')
    ap.add_argument('--top-n', type=int, default=15)
    ap.add_argument('--min-reads', type=int, default=5)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    rank_order = [r.strip().upper() for r in args.ranks.split(',') if r.strip()]
    rank_pos = {r: i for i, r in enumerate(rank_order)}

    # ---- flags for this sample ----
    call_of = {}
    if args.flags and os.path.exists(args.flags):
        with open(args.flags) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            try:
                si, ti, ci = header.index('sample'), header.index('taxid'), header.index('call')
            except ValueError:
                si = ti = ci = None
            if si is not None:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    if len(f) > max(si, ti, ci) and f[si] == args.sample:
                        call_of[f[ti]] = f[ci]

    nodes = parse_tree(args.report)

    def main_rank(code):
        l = norm_rank(code)
        return l if (l and l in rank_pos) else None

    # canonical candidate nodes (rank in requested set, >= min reads)
    canon = [i for i, n in enumerate(nodes)
             if main_rank(n['rank']) and n['clade_reads'] >= args.min_reads]

    # keep top-N per rank
    keep = set()
    by_rank = {}
    for i in canon:
        by_rank.setdefault(main_rank(nodes[i]['rank']), []).append(i)
    for r, idxs in by_rank.items():
        idxs.sort(key=lambda i: -nodes[i]['clade_reads'])
        keep.update(idxs[:args.top_n])

    def nearest_kept_ancestor(i):
        p = nodes[i]['parent']
        while p is not None:
            if p in keep:
                return p
            p = nodes[p]['parent']
        return None

    def domain_of(i):
        p = i
        seen = 0
        while p is not None and seen < 50:
            if norm_rank(nodes[p]['rank']) == 'D':
                return nodes[p]['name']
            p = nodes[p]['parent']
            seen += 1
        return 'unknown'

    # collapse to only the ranks actually present, so columns fill full width
    present = sorted({main_rank(nodes[i]['rank']) for i in keep}, key=lambda r: rank_pos[r])
    col_of = {r: k for k, r in enumerate(present)}
    ncol = max(len(present) - 1, 1)

    # ---- build plotly node/link arrays ----
    node_ids = {}         # keep-index -> plotly node index
    labels, node_colors, node_x = [], [], []
    for i in sorted(keep, key=lambda i: (rank_pos[main_rank(nodes[i]['rank'])], -nodes[i]['clade_reads'])):
        node_ids[i] = len(labels)
        n = nodes[i]
        rm = main_rank(n['rank'])
        labels.append(f"{n['name']} ({n['clade_reads']:,})")
        # colour: species by call if available, else by domain
        if rm == 'S' and n['taxid'] in call_of:
            node_colors.append(CALL_COLORS.get(call_of[n['taxid']], '#BFBFBF'))
        else:
            node_colors.append(DOMAIN_COLORS.get(domain_of(i), DOMAIN_COLORS['unknown']))
        node_x.append(0.01 + 0.98 * col_of[rm] / ncol)

    src, tgt, val, lcolor = [], [], [], []
    for i in keep:
        anc = nearest_kept_ancestor(i)
        if anc is None:
            continue
        src.append(node_ids[anc])
        tgt.append(node_ids[i])
        val.append(max(nodes[i]['clade_reads'], 1))
        base = DOMAIN_COLORS.get(domain_of(i), DOMAIN_COLORS['unknown'])
        lcolor.append(_rgba(base, 0.35))

    _render(args, labels, node_colors, node_x, src, tgt, val, lcolor, len(keep))
    print(f"[make_sankey] {args.sample}: {len(keep)} nodes, {len(src)} links",
          file=sys.stderr)


def _rgba(hexc, a):
    hexc = hexc.lstrip('#')
    r, g, b = int(hexc[0:2], 16), int(hexc[2:4], 16), int(hexc[4:6], 16)
    return f"rgba({r},{g},{b},{a})"


def _render(args, labels, node_colors, node_x, src, tgt, val, lcolor, nkeep):
    title = f"Taxonomic flow — {args.sample}"
    if nkeep == 0:
        with open(args.out, 'w') as fh:
            fh.write(f"<html><body style='font-family:sans-serif'>"
                     f"<h3>{title}</h3><p>No taxa above the reporting "
                     f"threshold (min {args.min_reads} clade reads).</p>"
                     f"</body></html>")
        return
    import plotly.graph_objects as go
    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node=dict(label=labels, color=node_colors, x=node_x,
                  pad=12, thickness=16,
                  line=dict(color='rgba(0,0,0,0.25)', width=0.5),
                  hovertemplate='%{label}<extra></extra>'),
        link=dict(source=src, target=tgt, value=val, color=lcolor,
                  hovertemplate='%{source.label} → %{target.label}<br>'
                                '%{value:,} reads<extra></extra>'),
    ))
    fig.update_layout(
        title=dict(text=title, x=0.02, font=dict(size=18)),
        font=dict(size=11), margin=dict(l=10, r=10, t=50, b=30), height=720,
        annotations=[dict(text=("Green = statistically supported · Grey = "
                                "likely contaminant · Flow width = reads"),
                          showarrow=False, x=0.02, y=-0.06, xref='paper',
                          yref='paper', font=dict(size=10, color='#666'))])
    fig.write_html(args.out, include_plotlyjs=True, full_html=True)


if __name__ == '__main__':
    main()
PYEOF

        python3 make_sankey.py \
            --report ~{report} \
            --sample ~{sample_id} \
            --flags ~{flags} \
            --ranks ~{ranks} \
            --top-n ~{top_n} \
            --min-reads ~{min_reads} \
            --out ~{sample_id}.sankey.html
    >>>

    output { File html = "~{sample_id}.sankey.html" }
    runtime {
        docker: "python:3.11-slim"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 20 HDD"
        preemptible: preemptible_tries
    }
}

task Summary {
    input {
        Array[File] sankeys
        Array[File] kronas
        File        significant
        File        summary_json
        Array[File] qc_jsons
        Array[File] host_stats
        Int         preemptible_tries
    }

    command <<<
        set -euo pipefail
        mkdir -p sankey krona qc host
        for f in ~{sep=' ' sankeys};    do cp "$f" sankey/ ; done
        for f in ~{sep=' ' kronas};     do cp "$f" krona/  ; done
        for f in ~{sep=' ' qc_jsons};   do cp "$f" qc/     ; done
        for f in ~{sep=' ' host_stats}; do cp "$f" host/   ; done

cat > make_summary.py <<'PYEOF'
#!/usr/bin/env python3
"""
Assemble a single index.html overview linking per-sample Sankey + Krona charts,
the significant-taxa table, QC stats (nanoq) and host-removal stats (flagstat).
Pure stdlib. Links are relative to results/summary/ (siblings: ../sankey, ../krona).
"""
import argparse
import glob
import json
import os
import re
import sys


def read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception:
        return {}


def parse_flagstat(path):
    total = mapped = 0
    try:
        with open(path) as fh:
            for line in fh:
                if 'in total' in line:
                    total = int(line.split()[0])
                elif re.search(r'\bmapped\b', line) and 'primary' not in line and '%' in line:
                    mapped = int(line.split()[0])
    except Exception:
        pass
    pct = (100.0 * mapped / total) if total else 0.0
    return total, mapped, pct


def sample_from(path, suffixes):
    b = os.path.basename(path)
    for s in suffixes:
        if b.endswith(s):
            return b[:-len(s)]
    return b.split('.')[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sankey-dir', default='sankey')
    ap.add_argument('--krona-dir', default='krona')
    ap.add_argument('--significant', required=True)
    ap.add_argument('--summary-json', required=True)
    ap.add_argument('--qc-dir', default='qc')
    ap.add_argument('--host-dir', default='host')
    ap.add_argument('--out', default='index.html')
    ap.add_argument('--table', default='run_summary.tsv')
    args = ap.parse_args()

    summary = read_json(args.summary_json)
    mode = summary.get('mode', 'n/a')

    sankeys = {sample_from(p, ['.sankey.html']): os.path.basename(p)
               for p in sorted(glob.glob(os.path.join(args.sankey_dir, '*.sankey.html')))}
    kronas = {sample_from(p, ['_krona.html']): os.path.basename(p)
              for p in sorted(glob.glob(os.path.join(args.krona_dir, '*_krona.html')))}
    qc = {sample_from(p, ['.qc.json']): read_json(p)
          for p in glob.glob(os.path.join(args.qc_dir, '*.qc.json'))}
    host = {sample_from(p, ['.host.flagstat']): parse_flagstat(p)
            for p in glob.glob(os.path.join(args.host_dir, '*.host.flagstat'))}

    samples = sorted(set(sankeys) | set(kronas) | set(qc) | set(host))

    # significant taxa table
    sig_rows = []
    try:
        with open(args.significant) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            for line in fh:
                sig_rows.append(dict(zip(header, line.rstrip('\n').split('\t'))))
    except Exception:
        header = []

    # ---- run_summary.tsv ----
    with open(args.table, 'w') as fh:
        fh.write('sample\tqc_reads\thost_pct\tn_significant\n')
        sig_by_sample = {}
        for r in sig_rows:
            sig_by_sample[r.get('sample', '')] = sig_by_sample.get(r.get('sample', ''), 0) + 1
        for s in samples:
            reads = qc.get(s, {}).get('reads', qc.get(s, {}).get('n', ''))
            hp = f"{host.get(s, (0, 0, 0))[2]:.2f}" if s in host else ''
            fh.write(f"{s}\t{reads}\t{hp}\t{sig_by_sample.get(s, 0)}\n")

    # ---- HTML ----
    css = """
    body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;color:#1c1c1e;background:#fafafa}
    header{background:#1f2d3d;color:#fff;padding:20px 28px}
    header h1{margin:0;font-size:20px} header p{margin:4px 0 0;opacity:.8;font-size:13px}
    .wrap{max-width:1100px;margin:0 auto;padding:22px}
    .banner{background:#fff3cd;border:1px solid #ffe69c;color:#664d03;padding:10px 14px;border-radius:8px;font-size:13px;margin-bottom:18px}
    .cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:14px}
    .card{background:#fff;border:1px solid #e5e5ea;border-radius:10px;padding:14px}
    .card h3{margin:0 0 8px;font-size:15px}
    .card a{display:inline-block;margin-right:10px;font-size:13px;color:#0a84ff;text-decoration:none}
    .metric{font-size:12px;color:#555;margin-top:6px}
    table{border-collapse:collapse;width:100%;background:#fff;font-size:13px;margin-top:8px}
    th,td{border:1px solid #e5e5ea;padding:6px 9px;text-align:left}
    th{background:#f2f2f7} tr:nth-child(even){background:#fafafe}
    .pill{padding:1px 7px;border-radius:10px;font-size:11px;color:#fff}
    .sig{background:#2ca02c} h2{font-size:16px;margin-top:26px}
    """
    h = [f"<html><head><meta charset='utf-8'><title>ont-metabuli summary</title>"
         f"<style>{css}</style></head><body>"]
    h.append("<header><h1>ont-metabuli — run summary</h1>"
             f"<p>Statistics mode: <b>{mode}</b> · "
             f"{summary.get('n_biological_samples','?')} samples · "
             f"{summary.get('n_controls',0)} negative control(s) · "
             f"{summary.get('n_significant_total','?')} significant taxa</p></header>")
    h.append("<div class='wrap'>")
    if summary.get('caveat'):
        h.append(f"<div class='banner'>⚠ {summary['caveat']}</div>")

    # per-sample cards
    h.append("<h2>Samples</h2><div class='cards'>")
    for s in samples:
        links = []
        if s in sankeys:
            links.append(f"<a href='../sankey/{sankeys[s]}'>Sankey ↗</a>")
        if s in kronas:
            links.append(f"<a href='../krona/{kronas[s]}'>Krona ↗</a>")
        metrics = []
        if s in host:
            metrics.append(f"host mapped: {host[s][2]:.1f}%")
        q = qc.get(s, {})
        if q:
            rd = q.get('reads', q.get('n'))
            if rd is not None:
                metrics.append(f"reads (QC): {rd}")
        nsig = sum(1 for r in sig_rows if r.get('sample') == s)
        metrics.append(f"significant taxa: {nsig}")
        h.append(f"<div class='card'><h3>{s}</h3>{' '.join(links) or '<i>no charts</i>'}"
                 f"<div class='metric'>{' · '.join(metrics)}</div></div>")
    h.append("</div>")

    # significant taxa table (top 50)
    h.append("<h2>Significant / detected taxa</h2>")
    if sig_rows:
        show = ['sample', 'name', 'reads', 'rel_abundance', 'fold_enrichment',
                'q_value', 'call']
        show = [c for c in show if c in header]
        h.append("<table><tr>" + ''.join(f"<th>{c}</th>" for c in show) + "</tr>")
        for r in sig_rows[:50]:
            cells = []
            for c in show:
                v = r.get(c, '')
                if c == 'call':
                    v = f"<span class='pill sig'>{v}</span>"
                cells.append(f"<td>{v}</td>")
            h.append("<tr>" + ''.join(cells) + "</tr>")
        h.append("</table>")
        if len(sig_rows) > 50:
            h.append(f"<p class='metric'>… {len(sig_rows)-50} more in "
                     f"stats/significant_taxa.tsv</p>")
    else:
        h.append("<p class='metric'>No taxa passed the significance thresholds.</p>")

    h.append("<h2>Outputs</h2><ul class='metric'>"
             "<li><code>metabuli/&lt;sample&gt;/</code> — raw classifications & report</li>"
             "<li><code>abundance/</code> — combined long table & count matrix</li>"
             "<li><code>stats/</code> — taxon_stats.tsv, significant_taxa.tsv, stats_summary.json</li>"
             "<li><code>sankey/</code> — interactive Sankey per sample</li>"
             "<li><code>pipeline_info/</code> — Nextflow timeline/report/trace</li></ul>")
    h.append("</div></body></html>")

    with open(args.out, 'w') as fh:
        fh.write('\n'.join(h))
    print(f"[make_summary] {len(samples)} samples, {len(sig_rows)} significant rows",
          file=sys.stderr)


if __name__ == '__main__':
    main()
PYEOF

        python3 make_summary.py \
            --sankey-dir sankey --krona-dir krona \
            --significant ~{significant} --summary-json ~{summary_json} \
            --qc-dir qc --host-dir host \
            --out index.html --table run_summary.tsv
    >>>

    output {
        File html  = "index.html"
        File table = "run_summary.tsv"
    }
    runtime {
        docker: "python:3.11-slim"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 HDD"
        preemptible: preemptible_tries
    }
}

task Flye {
    input {
        String sample_id
        String bin_label
        File   reads
        String read_type
        Int    min_bin_reads
        Int    preemptible_tries
        Int    cpu = 8
        Int    memory_gb = 32
    }
    Int disk_gb = ceil(size(reads, "GB") * 10) + 40

    command <<<
        set -euo pipefail
        out="~{sample_id}.~{bin_label}.assembly.fasta"
        n=$(( $(wc -l < ~{reads}) / 4 ))
        echo "[flye] ~{sample_id} ~{bin_label}: $n reads (minimum ~{min_bin_reads})"

        if [ "$n" -lt ~{min_bin_reads} ]; then
            echo "[flye] too few reads, skipping assembly"
            : > "$out"
            exit 0
        fi

        # --meta copes with the uneven coverage typical of a taxon bin; Flye
        # runs its own polishing iterations, so no separate polisher is needed.
        if flye --~{read_type} ~{reads} --meta --threads ~{cpu} --out-dir flye_out; then
            if [ -s flye_out/assembly.fasta ]; then
                cp flye_out/assembly.fasta "$out"
                [ -f flye_out/assembly_info.txt ] && cp flye_out/assembly_info.txt "~{sample_id}.~{bin_label}.assembly_info.txt" || true
            else
                : > "$out"
            fi
        else
            echo "[flye] assembly produced no contigs"
            : > "$out"
        fi
    >>>

    output { File assembly = "~{sample_id}.~{bin_label}.assembly.fasta" }
    runtime {
        docker: "quay.io/biocontainers/flye:2.9.6--py311h93bbee8_1"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible_tries
    }
}

task ShortReadAssembly {
    input {
        String      sample_id
        String      bin_label
        Array[File] reads          # [R1, R2]
        String      assembler      # megahit | spades
        Int         min_bin_reads
        Int         min_contig_len
        Int         preemptible_tries
        Int         cpu = 8
        Int         memory_gb = 32
    }
    Int disk_gb = ceil(size(reads, "GB") * 12) + 40

    command <<<
        set -euo pipefail
        out="~{sample_id}.~{bin_label}.assembly.fasta"
        n=$(( $(wc -l < ~{reads[0]}) / 4 ))
        echo "[~{assembler}] ~{sample_id} ~{bin_label}: $n pairs (minimum ~{min_bin_reads})"

        if [ "$n" -lt ~{min_bin_reads} ]; then
            echo "[~{assembler}] too few reads, skipping assembly"
            : > "$out"
            exit 0
        fi

        if [ "~{assembler}" = "spades" ]; then
            # metaSPAdes: higher contiguity, substantially more memory/time.
            # -m is a hard cap so the task cannot exceed its requested memory.
            if spades.py --meta -1 ~{reads[0]} -2 ~{reads[1]} \
                    -t ~{cpu} -m ~{memory_gb} -o asm_out; then
                if [ -s asm_out/contigs.fasta ]; then
                    awk -v m=~{min_contig_len} '/^>/{keep=0; n=$0; sub(/.*_length_/,"",n); sub(/_cov.*/,"",n); if (n+0>=m) keep=1} keep' \
                        asm_out/contigs.fasta > "$out" || cp asm_out/contigs.fasta "$out"
                else
                    : > "$out"
                fi
            else
                echo "[spades] assembly produced no contigs"; : > "$out"
            fi
        else
            # MEGAHIT: fast and memory-light. NOTE the bioconda osx-arm64 build
            # segfaults during assembly; the linux/amd64 image used here is fine.
            if megahit -1 ~{reads[0]} -2 ~{reads[1]} \
                    --min-contig-len ~{min_contig_len} -t ~{cpu} -o asm_out; then
                if [ -s asm_out/final.contigs.fa ]; then
                    cp asm_out/final.contigs.fa "$out"
                else
                    : > "$out"
                fi
            else
                echo "[megahit] assembly produced no contigs"; : > "$out"
            fi
        fi
    >>>

    output { File assembly = "~{sample_id}.~{bin_label}.assembly.fasta" }
    runtime {
        docker: if assembler == "spades"
                then "quay.io/biocontainers/spades:4.3.0--hde4eca7_1"
                else "quay.io/biocontainers/megahit:1.2.9--haf24da9_8"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible_tries
    }
}

task AssemblyStats {
    input {
        String sample_id
        String bin_label
        File   assembly
        Int    preemptible_tries
    }

    command <<<
        set -euo pipefail
        out="~{sample_id}.~{bin_label}.assembly_stats.tsv"
        if [ -s ~{assembly} ]; then
            seqkit stats -a -T ~{assembly} \
              | awk -v s='~{sample_id}' -v b='~{bin_label}' 'NR==1{print "sample\tbin\t"$0} NR>1{print s"\t"b"\t"$0}' > "$out"
        else
            printf 'sample\tbin\tnote\n%s\t%s\tno_assembly\n' '~{sample_id}' '~{bin_label}' > "$out"
        fi
    >>>

    output { File stats = "~{sample_id}.~{bin_label}.assembly_stats.tsv" }
    runtime {
        docker: "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 HDD"
        preemptible: preemptible_tries
    }
}

task ClassifyContigs {
    input {
        Array[File] assemblies
        File        db_tar
        Int         max_ram_gb
        Int         cpu
        Int         memory_gb
        Int         disk_gb
        String      disk_type
    }

    command <<<
        set -euo pipefail
        mkdir -p dbroot out reports
        tar -xzf ~{db_tar} -C dbroot
        DBDIR=$(dirname "$(find dbroot -maxdepth 4 -name db.parameters | head -1)")

        printf 'assembly\tclade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n' > reports/contig_reports.tsv
        for f in ~{sep=' ' assemblies}; do
            [ -s "$f" ] || continue
            b=$(basename "$f" .assembly.fasta)
            metabuli classify --seq-mode 3 --threads ~{cpu} --max-ram ~{max_ram_gb} \
                "$f" "$DBDIR" out "${b}_contigs" || continue
            if [ -f "out/${b}_contigs_report.tsv" ]; then
                grep -v '^#' "out/${b}_contigs_report.tsv" \
                  | awk -v b="$b" 'BEGIN{OFS="\t"} {print b, $0}' >> reports/contig_reports.tsv
            fi
        done
    >>>

    output { File merged_report = "reports/contig_reports.tsv" }
    runtime {
        docker: "quay.io/biocontainers/metabuli:1.2.0--pl5321h0bb26bb_0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " " + disk_type
        preemptible: 0
    }
}

task MapToReference {
    input {
        String        sample_id
        String        bin_label
        String        taxid
        Array[File]   reads        # [R1] or [R1, R2]
        Boolean       paired
        Array[String] ref_taxids
        Array[File]   references
        Int           cpu
        Int           memory_gb
        Int           preemptible_tries
    }
    Int disk_gb = ceil(size(reads, "GB") * 8 + size(references, "GB") * 6) + 30

    command <<<
        set -euo pipefail
        pfx="~{sample_id}.~{bin_label}"

        # pick the reference paired with this bin's taxid
        taxids=(~{sep=' ' ref_taxids})
        refs=(~{sep=' ' references})
        REF_SRC=""
        for i in "${!taxids[@]}"; do
            if [ "${taxids[$i]}" = "~{taxid}" ]; then REF_SRC="${refs[$i]}"; break; fi
        done
        if [ -z "$REF_SRC" ]; then
            echo "ERROR: no reference supplied for taxid ~{taxid}" >&2
            exit 1
        fi
        echo "reference for taxid ~{taxid}: $REF_SRC"

        # bcftools downstream needs a plain, faidx-able reference
        case "$REF_SRC" in
            *.gz) gunzip -c "$REF_SRC" > reference.fa ;;
            *)    cp "$REF_SRC" reference.fa ;;
        esac
        samtools faidx reference.fa

        minimap2 -ax ~{if paired then "sr" else "map-ont"} -t ~{cpu} reference.fa ~{sep=' ' reads} \
            | samtools sort -@ ~{cpu} -o "${pfx}.bam"
        samtools index "${pfx}.bam"

        samtools coverage "${pfx}.bam" > "${pfx}.coverage.txt"
        samtools consensus -f fasta -o "${pfx}.consensus.fasta" "${pfx}.bam" \
            || : > "${pfx}.consensus.fasta"
    >>>

    output {
        File bam                    = "~{sample_id}.~{bin_label}.bam"
        File bam_index              = "~{sample_id}.~{bin_label}.bam.bai"
        File consensus              = "~{sample_id}.~{bin_label}.consensus.fasta"
        File coverage               = "~{sample_id}.~{bin_label}.coverage.txt"
        File prepared_reference     = "reference.fa"
        File prepared_reference_fai = "reference.fa.fai"
    }
    runtime {
        docker: "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:b411340b52d82a9c276d87c7a3dcffc880be762f-0"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: preemptible_tries
    }
}

task CallVariants {
    input {
        String sample_id
        String bin_label
        File   bam
        File   bam_index
        File   reference
        File   reference_fai
        Int    preemptible_tries
        Int    cpu = 2
        Int    memory_gb = 8
    }
    Int disk_gb = ceil(size(bam, "GB") * 3 + size(reference, "GB") * 3) + 20

    command <<<
        set -euo pipefail
        pfx="~{sample_id}.~{bin_label}"
        # keep the reference and its index side by side for htslib
        cp ~{reference} reference.fa
        cp ~{reference_fai} reference.fa.fai

        bcftools mpileup -f reference.fa ~{bam} -Ou \
            | bcftools call -mv -Oz -o "${pfx}.variants.vcf.gz"
        bcftools index "${pfx}.variants.vcf.gz"
    >>>

    output {
        File vcf       = "~{sample_id}.~{bin_label}.variants.vcf.gz"
        File vcf_index = "~{sample_id}.~{bin_label}.variants.vcf.gz.csi"
    }
    runtime {
        docker: "quay.io/biocontainers/bcftools:1.24--h487d631_1"
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: preemptible_tries
    }
}
