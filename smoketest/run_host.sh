#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$JAVA_HOME/bin:$PWD/bin_ext/metabuli/bin:$HOME/miniforge3/envs/ont-metabuli/bin:$PATH"
export NXF_ANSI_LOG=false
echo "[$(date '+%H:%M:%S')] full pipeline WITH host removal (minimap2 vs CHM13)…"
nextflow run . -profile standard \
  --input smoketest/samplesheet_host.csv \
  --host_fasta references/chm13v2.0.fa.gz \
  --metabuli_db databases/refseq_standard/refseq_standard \
  --metabuli_max_ram 10 \
  --stats_min_reads 5 --stats_min_abundance 0.001 --sankey_min_reads 2 \
  --outdir smoketest/results_host -work-dir smoketest/work_host
echo "[$(date '+%H:%M:%S')] exit=$?"
echo "=== host removal stats ==="; cat smoketest/results_host/host_removed/HOSTTEST.host.flagstat 2>/dev/null | head -6
