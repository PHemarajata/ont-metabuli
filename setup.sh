#!/usr/bin/env bash
# ---------------------------------------------------------------------------
#  ont-metabuli setup — installs everything needed to run natively on macOS
#  (Apple Silicon or Intel) with minimal fuss.
#
#  Usage:
#    ./setup.sh                 # env + native Metabuli binary  (do this first)
#    ./setup.sh --host          #   + download T2T-CHM13v2 host reference
#    ./setup.sh --db RefSeq_virus   #   + download a prebuilt Metabuli DB
#    ./setup.sh --host --db RefSeq_virus   # everything
#
#  DB options (RAM-aware): RefSeq_virus (~8 GiB, fits 16 GB) | GTDB | RefSeq | RefSeq_release
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")"
ROOT="$(pwd)"
BIN="$ROOT/bin_ext"        # native Metabuli lives here
REF="$ROOT/references"
DBROOT="$ROOT/databases"
ENV_NAME="ont-metabuli"

DO_HOST=0; DB_NAME=""; DO_MEDAKA=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --host)   DO_HOST=1; shift;;
    --db)     DB_NAME="${2:-}"; shift 2;;
    --medaka) DO_MEDAKA=1; shift;;   # optional separate env for --polisher medaka
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

OS="$(uname -s)"; ARCH="$(uname -m)"
echo "==> Host: $OS / $ARCH"

# --- conda/mamba env ---------------------------------------------------------
CONDA_BIN="$(command -v mamba || command -v conda || true)"
if [[ -z "$CONDA_BIN" ]]; then
  echo "!! mamba/conda not found. Install miniforge: https://github.com/conda-forge/miniforge"
  exit 1
fi
echo "==> Creating/updating conda env '$ENV_NAME' (native)..."
if ! conda env list | grep -qE "^${ENV_NAME}\s"; then
  "$CONDA_BIN" env create -n "$ENV_NAME" -f env/environment.yml
else
  "$CONDA_BIN" env update -n "$ENV_NAME" -f env/environment.yml
fi

# --- native Metabuli binary --------------------------------------------------
# The new prebuilt databases (opendata.mmseqs.org) REQUIRE Metabuli >= v1.2.0.
# Pin to the versioned GitHub release asset (reproducible) rather than the
# rolling mmseqs.com build. v1.2.0 adds reduced-alphabet + spaced k-mers +
# syncmers, which the new DBs are built with.
METABULI_VERSION="1.2.0"
GH="https://github.com/steineggerlab/Metabuli/releases/download/${METABULI_VERSION}"
mkdir -p "$BIN"
if [[ "$OS" == "Darwin" ]]; then
  URL="$GH/metabuli-osx-universal.tar.gz"        # native on Apple Silicon + Intel
elif [[ "$OS" == "Linux" ]]; then
  ARCH_TAG="avx2"; case "$ARCH" in arm64|aarch64) ARCH_TAG="arm64";; esac
  URL="$GH/metabuli-linux-${ARCH_TAG}.tar.gz"
else
  echo "!! Unsupported OS for auto Metabuli install: $OS"; URL=""
fi
METABULI="$BIN/metabuli/bin/metabuli"
# (re)install if missing, or if the existing binary predates syncmer support (< v1.2.0)
need_install=0
if [[ ! -x "$METABULI" ]]; then need_install=1
elif ! "$METABULI" classify -h 2>&1 | grep -qi 'syncmer'; then
  echo "==> Existing Metabuli is older than v1.2.0 (no syncmer support); reinstalling."
  rm -rf "$BIN/metabuli"; need_install=1
fi
if [[ -n "$URL" && "$need_install" == "1" ]]; then
  echo "==> Downloading Metabuli v${METABULI_VERSION}: $URL"
  curl -fSL "$URL" -o "$BIN/metabuli.tar.gz"
  tar -xzf "$BIN/metabuli.tar.gz" -C "$BIN"
  rm -f "$BIN/metabuli.tar.gz"
fi
if [[ -x "$METABULI" ]]; then
  if "$METABULI" classify -h 2>&1 | grep -qi 'syncmer'; then
    echo "==> Metabuli ready (v${METABULI_VERSION}+, supports new databases)."
  else
    echo "!! WARNING: installed Metabuli lacks syncmer support — new DBs may fail."
  fi
fi

# --- host reference (T2T-CHM13v2) -------------------------------------------
# Contig names don't matter for host removal (mapped reads are discarded), so we
# prefer the NCBI mirror (fast + reliable). The human-pangenomics S3 "analysis
# set" is a fallback but has been throttling badly (~100 KB/s) and its old
# path-style URL 301-redirects to a 4 KB error file — use virtual-hosted style.
if [[ "$DO_HOST" == "1" ]]; then
  mkdir -p "$REF"
  HOST_FA="$REF/chm13v2.0.fa.gz"
  HOST_NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
  HOST_S3="https://human-pangenomics.s3.amazonaws.com/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
  # (re)download if missing or if a previous attempt left a tiny error page
  if [[ ! -f "$HOST_FA" ]] || [[ $(stat -f%z "$HOST_FA" 2>/dev/null || echo 0) -lt 1000000 ]]; then
    echo "==> Downloading T2T-CHM13v2 (~890 MB) from NCBI..."
    if ! curl -fSL -C - --retry 10 --retry-delay 5 --retry-connrefused "$HOST_NCBI" -o "$HOST_FA"; then
      echo "==> NCBI failed; trying human-pangenomics S3 (may be slow)..."
      curl -fSL -C - --retry 10 --retry-delay 5 "$HOST_S3" -o "$HOST_FA"
    fi
  fi
  # sanity: must be a valid gzip, not an XML/HTML error page
  if [[ "$(file -b "$HOST_FA" 2>/dev/null)" == gzip* ]] && gzip -t "$HOST_FA" 2>/dev/null; then
    echo "==> Host reference: $HOST_FA ($(du -h "$HOST_FA" | cut -f1))"
  else
    echo "!! Host download looks wrong (not a valid gzip): $HOST_FA"; head -c 200 "$HOST_FA"; echo
  fi
fi

# --- Metabuli database -------------------------------------------------------
# The binary's built-in `metabuli databases` command points at a decommissioned
# endpoint (steineggerlab.workers.dev -> S3 NoSuchKey). Databases now live at
# opendata.mmseqs.org. NOTE: as of 2026 there is NO small viral-only prebuilt DB
# anymore; the smallest curated DBs are tens of GB and run disk-backed via
# --metabuli_max_ram on a 16 GB machine. Valid names below.
DB_HOST="https://opendata.mmseqs.org/metabuli"
if [[ -n "$DB_NAME" ]]; then
  case "$DB_NAME" in
    gtdb226|gtdb232|refseq_standard|hrgm2|hrom) : ;;
    *) echo "!! Unknown DB '$DB_NAME'. Valid: hrom hrgm2 refseq_standard gtdb226 gtdb232"; exit 1;;
  esac
  mkdir -p "$DBROOT"
  TARBALL="$DBROOT/$DB_NAME.tar.gz"
  DB_URL="$DB_HOST/$DB_NAME.tar.gz"
  # resumable download: compare local vs remote size, resume (-C -) if partial
  REMOTE=$(curl -sIL "$DB_URL" | awk -F': ' 'tolower($1)=="content-length"{print $2}' | tr -d '\r' | tail -1)
  LOCAL=$(stat -f%z "$TARBALL" 2>/dev/null || echo 0)
  if [[ "$LOCAL" != "$REMOTE" ]]; then
    echo "==> Downloading Metabuli DB '$DB_NAME' from $DB_URL"
    echo "    (~$(( ${REMOTE:-0} / 1000000000 )) GB; resumable — re-run setup to continue if interrupted)"
    curl -fSL -C - --retry 8 --retry-delay 5 "$DB_URL" -o "$TARBALL"
  else
    echo "==> DB tarball already complete: $TARBALL"
  fi
  echo "==> Extracting $DB_NAME ..."
  mkdir -p "$DBROOT/$DB_NAME"
  tar -xzf "$TARBALL" -C "$DBROOT/$DB_NAME"
  # locate the directory that actually holds the Metabuli index files
  DBDIR=$(find "$DBROOT/$DB_NAME" -maxdepth 2 -type f \
            \( -name 'taxonomyDB' -o -name 'diffIdx' -o -name 'deltaIdx*' -o -name '*.mtbl' \) \
            2>/dev/null | head -1 | xargs -r dirname)
  DBDIR="${DBDIR:-$DBROOT/$DB_NAME}"
  echo "==> DB ready: $DBDIR"
  echo "    (pass this exact path to --metabuli_db; you may delete $TARBALL to reclaim space)"
fi

# --- optional medaka env (separate: python/samtools conflict on arm64) -------
if [[ "$DO_MEDAKA" == "1" ]]; then
  if ! conda env list 2>/dev/null | grep -qE "^ont-metabuli-medaka\s"; then
    echo "==> Creating separate 'ont-metabuli-medaka' env..."
    "$CONDA_BIN" create -n ont-metabuli-medaka -c bioconda -c conda-forge 'medaka>=2.1.1' -y || \
      echo "!! medaka env creation failed; use --polisher racon instead."
  fi
  MED_PREFIX="$(conda env list 2>/dev/null | awk '/ont-metabuli-medaka/{print $NF}')"
  echo "==> medaka env: ${MED_PREFIX:-<not created>}"
  echo "    run with: --polisher medaka --medaka_env ${MED_PREFIX:-PATH} --medaka_model <model>"
fi

# --- write env.sh helper -----------------------------------------------------
cat > "$ROOT/env.sh" <<EOF
# source this before running the pipeline:  source env.sh
export PATH="$BIN/metabuli/bin:\$PATH"
# activate the native tools env (nanoq, minimap2, samtools, python libs):
#   conda activate $ENV_NAME
EOF

cat <<EOF

============================================================
 Setup complete.

 Next:
   conda activate $ENV_NAME
   source env.sh                 # puts native metabuli on PATH

 Then run (see README):
   nextflow run . -profile standard \\
     --input samplesheet.csv \\
     --host_fasta references/chm13v2.0.fa.gz \\
     --metabuli_db databases/${DB_NAME:-RefSeq_virus} \\
     --outdir results
============================================================
EOF
