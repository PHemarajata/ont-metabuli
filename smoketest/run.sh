#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$JAVA_HOME/bin:$PWD/bin_ext/metabuli/bin:$HOME/miniforge3/envs/ont-metabuli/bin:$PATH"
export NXF_ANSI_LOG=false
echo "[$(date '+%H:%M:%S')] tools: metabuli=$(command -v metabuli) chopper=$(command -v chopper) python3=$(command -v python3)"
echo "[$(date '+%H:%M:%S')] starting nextflow (real DB, disk-backed metabuli)…"
nextflow run . -profile standard \
  --input smoketest/samplesheet.csv \
  --metabuli_db databases/refseq_standard/refseq_standard \
  --skip_host_removal \
  --metabuli_max_ram 10 \
  --stats_min_reads 5 --stats_min_abundance 0.001 --sankey_min_reads 2 \
  --outdir smoketest/results \
  -work-dir smoketest/work
echo "[$(date '+%H:%M:%S')] nextflow exit=$?"
