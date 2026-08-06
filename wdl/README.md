# ont-metabuli — Terra / WDL

A WDL port of the [ont-metabuli Nextflow pipeline](../README.md) for
[Terra.bio](https://terra.bio). Same science, re-architected for Cromwell on GCP.

```
 sample_set
     │
     ▼
 ┌─────────┐   ┌──────────────┐   ┌──────────────────────────┐   ┌──────────────┐
 │ chopper │──▶│  minimap2    │──▶│  Metabuli classify       │──▶│ combine +    │
 │ QC      │   │ host removal │   │  + extract (ONE task,    │   │ abundance    │
 └─────────┘   └──────────────┘   │    one DB localization)  │   └──────┬───────┘
   scattered      scattered       └────────────┬─────────────┘          │
                                               │                 ┌──────▼───────┐
                                   optional    │                 │ decontam     │
                                   assembly ◀──┘                 │ statistics   │
                                   (Flye / consensus)            └──────┬───────┘
                                                          ┌─────────────┼──────────┐
                                                          ▼             ▼          ▼
                                                      Sankey         Krona    summary.html
```

## Files

| File | Purpose |
|------|---------|
| `ont_metabuli.wdl` | **The workflow.** Single self-contained file — import this into Terra. |
| `stage_database.wdl` | Run **once** to pull a Metabuli database into your workspace bucket. |
| `inputs.example.json` | Curated, commented starting point. |
| `inputs.template.json` | Full auto-generated input list (`womtool inputs`). |
| `ont_metabuli.wdl.tmpl` + `build_wdl.py` | Source template and generator — see [Regenerating](#regenerating-the-wdl). |
| `test/` | Local Cromwell smoke test with the database tasks stubbed. |

---

## 1. Stage the database (once)

The prebuilt databases live at `https://opendata.mmseqs.org/metabuli` and **all
require Metabuli ≥ 1.2.0** (the container pinned here is 1.2.0). Rather than
downloading ~74 GB to a laptop and pushing it back up, run the helper workflow
so the data never leaves Google's network:

1. Import `stage_database.wdl` into your workspace.
2. Run it with `db_name = "refseq_standard"` (no other inputs needed).
3. Copy the `database_tar` output path — that `gs://…` URI is your
   `ont_metabuli.metabuli_db_tar` from now on.

| `db_name` | Download | Contents |
|-----------|----------|----------|
| `refseq_standard` | ~75 GB | archaea, bacteria, virus, plasmid, protozoa, fungi + human |
| `hrom` | ~22 GB | human **oral** microbiome + human + virus |
| `hrgm2` | ~85 GB | human **gut** microbiome v2 + human + virus |
| `gtdb226` / `gtdb232` | ~378 GB | full GTDB + human + virus |

> There is no small viral-only database any more — the former `RefSeq_virus`
> (~4 GB) has been withdrawn from every host.

You will also want a host genome in the bucket, e.g. T2T-CHM13v2:

```bash
gsutil cp chm13v2.0.fa.gz gs://YOUR-WORKSPACE-BUCKET/reference/
```

---

## 1b. Platform: ONT or Illumina

Set `platform` per **run** — not per sample. The statistics compare samples
against each other, and ONT vs Illumina differ enough in bias and depth that
mixing them in one significance test would be confounded; run each separately.

| Stage | `platform = "ont"` (default) | `platform = "illumina"` |
|-------|------------------------------|--------------------------|
| QC | chopper | fastp (paired-end) |
| Host removal | minimap2 `map-ont` | minimap2 `sr`, keeps pairs with **both** mates unmapped |
| Classification | Metabuli `--seq-mode 3` | Metabuli `--seq-mode 2 --precise 1` |
| Assembly (optional) | Flye `--meta` | MEGAHIT (default) or metaSPAdes (`assembler = "spades"`) |

For Illumina, supply `fastqs_2` alongside `fastqs` (same order):

```
ont_metabuli.platform  = "illumina"
ont_metabuli.fastqs    = this.samples.fastq_1
ont_metabuli.fastqs_2  = this.samples.fastq_2
```

See `inputs.illumina.example.json`. Everything downstream — abundance tables,
statistics, Sankey, summary — is identical for both platforms.

## 2. Set up the data table

Create a `sample` table with at least `sample_id` and a FASTQ column, then group
the samples into a `sample_set`:

| sample_id | fastq | sample_type | batch |
|-----------|-------|-------------|-------|
| SAMPLE01 | gs://…/SAMPLE01.fastq.gz | sample | run1 |
| SAMPLE02 | gs://…/SAMPLE02.fastq.gz | sample | run1 |
| NTC01 | gs://…/NTC01.fastq.gz | negative_control | run1 |

Run the workflow **on the sample_set** (not per sample) — the statistics compare
samples against each other and against controls, so they need the whole batch at
once. Wire the inputs up as:

```
ont_metabuli.sample_ids   = this.samples.sample_id
ont_metabuli.fastqs       = this.samples.fastq
ont_metabuli.sample_types = this.samples.sample_type
ont_metabuli.batches      = this.samples.batch
```

`sample_types` and `batches` are optional. Any of
`control / negative_control / ntc / blank / neg` marks a negative control, which
switches the statistics into contaminant-aware mode (see
[the main README](../README.md#4-the-significance-statistics)).

---

## 3. Run

Import `ont_metabuli.wdl`, point `metabuli_db_tar` and `host_fasta` at the files
you staged, and launch. Start from `inputs.example.json`.

Key outputs: `summary_html` (start here), `sankey_html`, `significant_taxa`,
`taxon_stats`, `abundance_matrix`, `metabuli_reports`.

### Optional assembly on binned reads

Off by default; see [the main README](../README.md#4b-optional-assembly-on-binned-reads)
for the rationale.

| Input | Mode |
|-------|------|
| `assemble_unclassified = true` | de novo Flye on **unclassified** reads — the novelty hunt |
| `assemble_taxa = ["562","10239"]` | de novo draft genome per taxon |
| `consensus_taxids` + `consensus_references` | reference consensus + variants (positional pairs) |
| `classify_contigs = true` | re-classify contigs to name them (**re-localizes the database** — see costs) |

---

## Resources and cost

The database dominates everything. Localizing ~74 GB and unpacking it to ~109 GB
is the single most expensive step, which is why **classification and read
binning happen in one batched task** rather than scattered per sample.

| Task | Default shape | Notes |
|------|---------------|-------|
| `MetabuliClassifyAndExtract` | 16 CPU / 120 GB RAM / 300 GB SSD, non-preemptible | one DB localization for the whole run |
| `Minimap2Index` | 8 CPU / 32 GB / SSD | once per run, skipped if you pass `host_mmi` |
| `HostRemoval` | 8 CPU / 32 GB | scattered per sample |
| `Chopper`, Python tasks | 1–4 CPU / 4–8 GB | cheap, preemptible |

Practical guidance:

- **Keep `metabuli_max_ram_gb` below `classify_memory_gb`** (defaults: 100 vs 120).
  With the index resident in RAM this runs far faster than the disk-backed mode
  you get on a 16 GB laptop.
- **Adding samples is cheap; adding runs is not.** The DB localization is paid
  once per workflow submission, so batch your samples into one sample_set.
- **`classify_contigs` pays the DB cost a second time.** Only enable it when you
  actually need contigs named (mainly novelty mode).
- `classify_disk_gb` must exceed *tarball + extracted* (74 + 109 ≈ 185 GB for
  `refseq_standard`); the 300 GB default leaves headroom. Bump it for GTDB.
- The classify task is **non-preemptible** on purpose — losing it late would
  repeat the whole localization.

---

## Differences from the Nextflow pipeline

These are deliberate adaptations, not omissions:

| Nextflow | Terra/WDL | Why |
|----------|-----------|-----|
| classification scattered per sample | one batched task | avoids N × 74 GB DB localization |
| `--assemble_taxa significant` | **not ported** — list explicit taxids | binning needs the DB, so it runs inside the classify task, before statistics exist |
| `--polisher racon\|medaka` | Flye's built-in polishing | avoids a second polishing container; Flye already polishes. Medaka can be added if you need ONT neural consensus. |
| `abricate` AMR screen | not included | it was excluded natively for an arm64 perl dependency; on Terra run it separately over `assembly_fasta` |
| `-resume` | Cromwell call caching | enable "Use call caching" in the Terra launch dialog |
| samplesheet CSV | Terra data table / parallel arrays | idiomatic for Terra |

Everything else — QC thresholds, host-removal logic, Metabuli invocation
(`--seq-mode 3`), the statistics (Fisher/Poisson + BH), the Sankey and the HTML
summary — is identical, and the Python is byte-identical (see below).

---

## Regenerating the WDL

`ont_metabuli.wdl` is **generated**. The analysis scripts are embedded verbatim
so Terra needs no auxiliary file staging, and generating them guarantees the
embedded copies match the ones the Nextflow pipeline validated:

```bash
python3 wdl/build_wdl.py     # ont_metabuli.wdl.tmpl + bin/*.py -> ont_metabuli.wdl
```

Edit `ont_metabuli.wdl.tmpl` (workflow logic) or `../bin/*.py` (analysis), never
the generated file.

---

## Testing locally

```bash
java -jar womtool.jar validate wdl/ont_metabuli.wdl

python3 wdl/test/make_stub_wdl.py          # stubs the two DB-dependent tasks
cd wdl/test && java -jar cromwell.jar run ont_metabuli.stub.wdl -i inputs.local.json
```

The stub replaces only `MetabuliClassifyAndExtract` and `ClassifyContigs` with
fixtures that emit realistic Metabuli v1.2.0 output (full rank names, `#` header,
sub-ranks). Everything else — scatter routing, glob collection, bin labelling,
file staging and the real analysis scripts — executes for real. It is the WDL
counterpart of `nextflow -stub-run`.

## Verification status

What has actually been checked, and what has not:

- `womtool validate` passes for `ont_metabuli.wdl`, `stage_database.wdl`, and
  `inputs.example.json` (no unexpected-input warnings).
- **Every container image was pulled and its tools run**, including confirmation
  that `metabuli:1.2.0--pl5321h0bb26bb_0` reports version 1.2.0 and exposes
  `--syncmer` (required by the current databases) plus the `extract` options used
  here, and that the mulled image really provides minimap2 2.31 + samtools 1.23.1
  with the `consensus` and `coverage` subcommands.
- The embedded Python was executed inside `python:3.11-slim` with the **pinned**
  dependency versions against real Metabuli reports, reproducing the Nextflow
  results (the pins differ from the versions used natively, so this mattered).
- Metabuli's **paired-end behaviour was verified against the real database**, not
  assumed: `--seq-mode 2` with two query files classifies a pair as one query,
  and `extract` emits two files named `<R1base>_<taxid>.fq` / `<R2base>_<taxid>.fq`
  (which is why the task stages inputs under canonical names).
- The stubbed workflow ran to `Succeeded` under Cromwell 92 locally in three
  configurations:
  - **full (ONT)** — host removal plus all three assembly modes: 3 samples
    produced 3 reports/Sankeys/Krona/QC/flagstats, 6 de novo assemblies (3 × 2
    bins) and 3 consensus + VCF sets, routed correctly (Flye vs consensus);
  - **Illumina** — same shape, with `Fastp` and `ShortReadAssembly` running in
    place of `Chopper` and `Flye`, and paired bins pairing correctly by index;
  - **defaults** — every optional module off and host removal skipped: core
    outputs populated and all disabled paths returning empty arrays rather than
    failing on empty globs/scatters.

Not verified: **a real Terra run against a real database.** Runtime and cost
figures above are estimates from the file sizes, not measurements, and the
Metabuli/ClassifyContigs command lines were exercised against the real container
CLI but not against a real index.
