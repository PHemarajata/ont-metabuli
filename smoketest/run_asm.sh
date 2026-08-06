#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$JAVA_HOME/bin:$PWD/bin_ext/metabuli/bin:$HOME/miniforge3/envs/ont-metabuli/bin:$PATH"
export NXF_ANSI_LOG=false
echo "[$(date '+%H:%M:%S')] tools: flye=$(command -v flye) racon=$(command -v racon) seqkit=$(command -v seqkit)"
echo "[$(date '+%H:%M:%S')] full pipeline + assembly (de novo taxon 561 + ref consensus vs lambda + unclassified)…"
nextflow run . -profile standard \
  --input smoketest/samplesheet_asm.csv \
  --metabuli_db databases/refseq_standard/refseq_standard \
  --skip_host_removal --metabuli_max_ram 10 \
  --assemble_taxa 561 \
  --reference_consensus "561:$PWD/smoketest/lambda.fa" \
  --assemble_unclassified \
  --polisher racon --min_bin_reads 20 \
  --stats_min_reads 5 --stats_min_abundance 0.001 --sankey_min_reads 2 \
  --outdir smoketest/results_asm -work-dir smoketest/work_asm
echo "[$(date '+%H:%M:%S')] exit=$?"
