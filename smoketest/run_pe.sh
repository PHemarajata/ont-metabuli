#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$JAVA_HOME/bin:$PWD/bin_ext/metabuli/bin:$HOME/miniforge3/envs/ont-metabuli/bin:$PATH"
export NXF_ANSI_LOG=false
echo "[$(date '+%H:%M:%S')] tools: fastp=$(command -v fastp) megahit=$(command -v megahit)"
echo "[$(date '+%H:%M:%S')] REAL Illumina paired-end run vs refseq_standard…"
nextflow run . -profile standard \
  --platform illumina \
  --input smoketest/samplesheet_pe.csv \
  --metabuli_db databases/refseq_standard/refseq_standard \
  --skip_host_removal --metabuli_max_ram 10 \
  --assemble_taxa 561 --assemble_unclassified --assembler spades \
  --reference_consensus "561:$PWD/smoketest/lambda.fa" \
  --min_bin_reads 20 --stats_min_reads 5 --stats_min_abundance 0.001 --sankey_min_reads 2 \
  --outdir smoketest/results_pe -work-dir smoketest/work_pe
echo "[$(date '+%H:%M:%S')] exit=$?"
