# ont-metabuli

A **low-resource Nextflow pipeline for Oxford Nanopore (ONT) metagenomics** built
around the [Metabuli](https://github.com/steineggerlab/Metabuli) classifier.
Designed to run on a laptop (tested target: **Apple Silicon Mac, 16 GB RAM**).

```
 samplesheet.csv
      │
      ▼
 ┌──────────┐   ┌──────────────┐   ┌──────────────┐   ┌───────────────┐
 │ chopper  │──▶│  minimap2    │──▶│   Metabuli   │──▶│  combine +    │
 │ QC/filter│   │ host removal │   │  classify    │   │  abundance    │
 └──────────┘   │ (map-ont,    │   │ (--seq-mode3)│   └──────┬────────┘
                │  keep unmap) │   └──────────────┘          │
                └──────────────┘                     ┌───────▼────────┐
                                                     │ decontam stats │
                                                     │ (NC or no-NC)  │
                                                     └───────┬────────┘
                                              ┌──────────────┼─────────────┐
                                              ▼              ▼             ▼
                                          Sankey (html)   Krona        summary/index.html
```

## Platforms

Set `--platform` per **run** (not per sample — comparing ONT against Illumina in
one significance test would be confounded by their different biases and depths):

| Stage | `--platform ont` (default) | `--platform illumina` |
|-------|---------------------------|------------------------|
| QC | **chopper** (length + mean quality) | **fastp** (paired-end adapter trim + quality) |
| Host removal | minimap2 `map-ont`, keep unmapped | minimap2 `sr`, keep pairs with **both** mates unmapped |
| Classification | Metabuli `--seq-mode 3` | Metabuli `--seq-mode 2 --precise 1` (short-read preset) |
| Assembly (optional) | **Flye** `--meta` | **MEGAHIT** (default) or **metaSPAdes** (`--assembler spades`) |

Everything downstream — abundance tables, significance statistics, Sankey and the
HTML summary — is identical for both.

## Why these tool choices (efficiency-first)

| Step | Tool | Why it's light |
|------|------|----------------|
| QC / filtering | **chopper** / **fastp** | Rust/C++, streaming, ~zero RAM (both native arm64) |
| Host removal | **minimap2** + samtools | preset-matched (`map-ont`/`sr`); memory bounded by `-I`; keeps only *unmapped* reads |
| Classification | **Metabuli** (native macOS universal binary) | runs natively on M-series (no Rosetta); `--max-ram` caps memory |
| Stats | Python + scipy | Fisher / Poisson tests, BH-FDR (no heavy deps) |
| Sankey | Python + plotly | self-contained interactive HTML (works offline) |

> **MEGAHIT on Apple Silicon:** the bioconda `osx-arm64` build **segfaults**
> during assembly (verified: exit −11 on any input) and silently yields empty
> contigs. The pipeline warns you if you hit this. Locally use
> `--assembler spades` (verified working, assembled a complete lambda genome) or
> run with `-profile docker`. The linux/amd64 container — including on Terra — is
> unaffected.

---

## 1. Install

```bash
cd ont_metagenomics
./setup.sh --host --db hrom       # env + native Metabuli v1.2.0 + CHM13 host + a DB
conda activate ont-metabuli
source env.sh                     # puts the native metabuli binary on PATH
```

`setup.sh` will:
- create a native conda env `ont-metabuli` (chopper, minimap2, samtools, python libs — all arm64),
- download the pinned **native Metabuli v1.2.0 universal binary** (no emulation),
- optionally download the **T2T-CHM13v2** host genome (`--host`),
- optionally download a **prebuilt Metabuli DB** (`--db NAME`) and print the exact
  DB path to pass to `--metabuli_db`.

### Choosing a database (the landscape changed in 2026)

Metabuli's prebuilt DBs now require **v1.2.0+** and are hosted at
`https://opendata.mmseqs.org/metabuli`. **There is no small viral-only DB anymore.**
The remaining curated DBs are large; on 16 GB RAM they run *disk-backed* via
`--metabuli_max_ram` (correct, just slower). Download sizes (compressed):

| DB (`--db`) | Download | Contents | 16 GB RAM |
|-------------|----------|----------|-----------|
| `hrom` | ~22 GB | Human **oral** microbiome + human + virus | disk-backed (smallest) |
| `refseq_standard` | ~75 GB | RefSeq archaea/bacteria/virus/fungi/protozoa + human | disk-backed |
| `hrgm2` | ~85 GB | Human **gut** microbiome + human + virus | disk-backed |
| `gtdb226` / `gtdb232` | ~378 GB | full GTDB (+ human + virus) | ❌ not practical |

All bundle the human T2T-CHM13v2 genome, but this pipeline still does a **dedicated
upfront host-removal** pass (minimap2) so host reads never reach the classifier.

> **Want fast, RAM-resident runs?** Build a **small custom DB** (e.g. viral or a
> targeted pathogen set) with `metabuli build` — a few GB, fits in RAM, no
> disk-backing. Ask and this repo can include a build helper.

> **Metabuli on ARM:** the bioconda `metabuli` package is x86-only. `setup.sh`
> avoids that by pinning the **native macOS universal binary** from the GitHub
> v1.2.0 release. Keep it on `PATH` (via `source env.sh`) and run `-profile standard`.

---

## 2. Prepare a samplesheet

**ONT** — one FASTQ per sample (see `assets/samplesheet_example.csv`):

```csv
sample,fastq,sample_type,batch
SAMPLE01,/abs/path/SAMPLE01.fastq.gz,sample,run1
SAMPLE02,/abs/path/SAMPLE02.fastq.gz,sample,run1
NTC01,/abs/path/NTC01.fastq.gz,negative_control,run1
```

**Illumina paired-end** — add a `fastq_2` column and pass `--platform illumina`
(see `assets/samplesheet_illumina_example.csv`):

```csv
sample,fastq,fastq_2,sample_type,batch
SAMPLE01,/abs/path/S1_R1.fastq.gz,/abs/path/S1_R2.fastq.gz,sample,run1
NTC01,/abs/path/NTC_R1.fastq.gz,/abs/path/NTC_R2.fastq.gz,negative_control,run1
```

- **`sample_type`** — `sample` for biological samples; any of
  `control,negative_control,ntc,blank,neg` marks a negative control.
  Controls are **optional** (leave them out if you have none).
- **`batch`** — optional; controls are paired to samples of the same batch
  (falls back to pooling all controls). Leave blank for a single pool.

---

## 3. Run

```bash
nextflow run . -profile standard \
    --input samplesheet.csv \
    --host_fasta references/chm13v2.0.fa.gz \
    --metabuli_db databases/hrom \
    --metabuli_max_ram 10 \
    --outdir results

# use the exact DB path setup.sh printed after "==> DB ready:"
```

Open **`results/summary/index.html`** when it finishes.

For **Illumina** paired-end data:

```bash
nextflow run . -profile standard \
    --platform illumina \
    --input samplesheet_illumina.csv \
    --host_fasta references/chm13v2.0.fa.gz \
    --metabuli_db databases/refseq_standard/refseq_standard \
    --assembler spades \
    --outdir results
```

Handy flags:
- `--platform ont|illumina` — selects the whole tool chain (see table above)
- `--min_qual 10 --min_len 500` — ONT QC thresholds
- `--sr_min_qual 15 --sr_min_len 50` — Illumina QC thresholds
- `--stats_rank species|genus`
- `--skip_host_removal` / `--skip_qc`
- `--max_cpus 8 --max_memory 12.GB` — resource ceilings
- `-resume` — continue from the last checkpoint

---

## 4. The significance statistics

Auto-selected from your samplesheet:

### With a negative control (rigorous)
For each taxon, its proportion of classified reads in the **sample** is compared to
the **batch-matched control** via a one-sided **Fisher exact test** (enriched in
sample), with **Benjamini–Hochberg** FDR across taxa. Calls:
- `contaminant` — proportion in control ≥ in sample
- `significant` — `q < alpha` **and** `fold ≥ min_fold`
- `below_threshold` / `not_significant`

Tunables: `--stats_alpha 0.05 --stats_min_fold 5 --stats_min_reads 10 --stats_min_abundance 1e-4`.

### Without a negative control (heuristic, honest)
There is no empirical blank to subtract, so significance is limited. The pipeline
reports a transparent detection heuristic: a one-sided **Poisson tail test** of each
taxon's read count against an assumed low background rate
(`--stats_bg_rate × classified_reads`), BH-corrected, plus a **Wilson 95% CI** on
relative abundance. Calls: `detected` / `below_threshold` / `not_significant`.

> These are **detection-confidence flags, not contamination-corrected results.**
> A blank/NTC is strongly recommended for any claim of "real" presence.

Outputs in `results/stats/`: `taxon_stats.tsv`, `significant_taxa.tsv`,
`stats_summary.json`.

---

## 4b. Optional: assembly on binned reads

Off by default. Metabuli's `extract` bins reads by taxon (or pulls the
*unclassified* pile), and this module assembles those bins — far lighter than
whole-metagenome co-assembly and aimed at a laptop. Three independent modes:

| Flag | Mode | Output |
|------|------|--------|
| `--assemble_unclassified` | **de novo on unclassified reads** (find novel/divergent agents) | `assembly/<sample>/unclassified/…polished.fasta` |
| `--assemble_taxa '562,10239'` or `'significant'` | **de novo per taxon** (draft genome of a hit) | `assembly/<sample>/taxon_<id>/…polished.fasta` |
| `--reference_consensus '562:/path/ref.fa;…'` | **reference consensus + variants** for a named organism | `reference_consensus/<sample>/ref_<id>/…consensus.fasta`, `…variants.vcf.gz` |

```bash
nextflow run . -profile standard --input samplesheet.csv \
    --metabuli_db databases/refseq_standard/refseq_standard \
    --host_fasta references/chm13v2.0.fa.gz \
    --assemble_taxa significant \
    --assemble_unclassified \
    --reference_consensus '10760:/refs/phage_lambda.fa' \
    --outdir results
```

- **Assembler:** chosen from `--platform`, override with `--assembler`:
  - ONT → **metaFlye** (`--flye_read_type nano-hq|nano-raw|nano-corr`)
  - Illumina → **MEGAHIT** (default) or **metaSPAdes** (`--assembler spades`,
    `--min_contig_len 500`)

  Bins under `--min_bin_reads` (default 50) are skipped — an unassemblable bin is
  normal, not a failure. Mixing a platform with the wrong assembler is rejected
  up front with a clear message.
- **Polishing:** `--polisher racon` (default; fast, no model) · `none` · `medaka`.
  **medaka** needs a separate env (python/samtools conflict on arm64): create it
  with `./setup.sh --medaka`, then run `--polisher medaka --medaka_env <prefix> --medaka_model <model>`.
- `--classify_contigs` re-runs Metabuli on the contigs to name them (rescans the
  DB — mainly useful to identify novelty-mode contigs).
- **AMR screening:** `abricate` has an unmet perl dependency on Apple Silicon.
  Run AMR/virulence on the produced `…polished.fasta` via `-profile docker`, on
  Linux, or a separate tool — it's intentionally not wired into the native path.

Assembly stays lightweight because each bin is a small read subset; a per-taxon
draft of one organism assembles in seconds–minutes on this hardware.

## 5. Outputs

```
results/
├── qc/                      chopper QC JSON reports
├── host_removed/            *.host.flagstat (host-mapping rate)
├── metabuli/<sample>/       *_report.tsv, *_classifications.tsv, *_krona.html
├── abundance/               combined_abundance_long.tsv, abundance_matrix_*.tsv
├── stats/                   taxon_stats.tsv, significant_taxa.tsv, stats_summary.json
├── sankey/                  <sample>.sankey.html   (interactive)
├── summary/index.html       ← start here
└── pipeline_info/           timeline / report / trace / DAG
```

Sankey colouring: **green** = statistically supported, **grey** = likely
contaminant, flow width = clade reads.

---

## 6. Validate the wiring (no data / DB needed)

```bash
nextflow run . -profile test -stub-run --outdir test_out
```

This runs every process with tiny stub commands to confirm the DAG connects
end-to-end. For a real smoke test, point `--metabuli_db` / `--host_fasta` at
small real paths.

---

## Profiles

- **`standard`** (recommended on Mac) — tools from PATH (native env + native Metabuli).
- **`conda`** — Nextflow builds a conda env per process (portable; Metabuli x86 only).
- **`docker`** — containers; on Apple Silicon runs under `linux/amd64` emulation.

## Citations

Please cite **Metabuli** (Kim et al.), **minimap2** (Li 2018), **samtools**,
**chopper**, **Nextflow** (Di Tommaso et al.), and **plotly** when publishing.

## License

Pipeline code: MIT. Bundled tools retain their own licenses.
