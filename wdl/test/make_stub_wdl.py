#!/usr/bin/env python3
"""
Produce a stubbed copy of ont_metabuli.wdl for local Cromwell testing.

The two tasks that need the ~74 GB Metabuli database (MetabuliClassifyAndExtract
and ClassifyContigs) have their command bodies replaced with fixtures that emit
realistically-formatted Metabuli v1.2.0 output. Everything else — scatter
routing, glob collection, sub()-based bin labelling, file staging and the real
analysis scripts — runs exactly as it would on Terra.

This is the WDL counterpart of `nextflow -stub-run`: it validates the plumbing,
not the classifier (whose CLI is verified separately against the real container).

Usage:  python3 wdl/test/make_stub_wdl.py  ->  wdl/test/ont_metabuli.stub.wdl
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WDL = os.path.join(os.path.dirname(HERE), 'ont_metabuli.wdl')
OUT = os.path.join(HERE, 'ont_metabuli.stub.wdl')

# A miniature Metabuli v1.2.0 report: full rank names, '#'-prefixed header,
# two-space indentation per level, and sub-ranks that must be ignored.
STUB_REPORT = r'''
cat > "reports/${sid}_report.tsv" <<'REPEOF'
#clade_proportion	clade_count	taxon_count	rank	taxID	name
100.0000	1000	20	no rank	1	root
80.0000	800	0	cellular root	131567	  cellular organisms
80.0000	800	5	domain	2	    Bacteria
60.0000	600	3	phylum	1224	      Pseudomonadota
55.0000	550	2	class	1236	        Gammaproteobacteria
50.0000	500	4	order	91347	          Enterobacterales
48.0000	480	3	family	543	            Enterobacteriaceae
45.0000	450	5	genus	561	              Escherichia
40.0000	400	400	species	562	                Escherichia coli
3.0000	30	30	strain	83333	                  Escherichia coli K-12
15.0000	150	150	species	1280	              Staphylococcus aureus
10.0000	100	0	domain	10239	    Viruses
9.0000	90	90	species	11676	      Human immunodeficiency virus 1
REPEOF
'''

STUB_CLASSIFY = r'''
        set -euo pipefail
        echo "[STUB] MetabuliClassifyAndExtract — database not localized"
        mkdir -p reports bins
        ids=(~{sep=' ' sample_ids})
        rds=(~{sep=' ' reads})
        PAIRED=~{if paired then "1" else "0"}

        for i in "${!ids[@]}"; do
            sid="${ids[$i]}"
            echo "[STUB] classifying $sid (paired=$PAIRED)"
''' + STUB_REPORT + r'''
            printf '#is_classified\tname\ttaxID\tquery_length\tscore\te_value\trank\n1\tr1\t562\t1500\t0.9\t0\tspecies\n' \
                > "reports/${sid}_classifications.tsv"
            echo "<html>krona ${sid}</html>" > "reports/${sid}_krona.html"

            emit_bin () {   # $1 label ; bins mirror the real naming convention
                printf '@r1\nACGTACGTACGTACGTACGTACGTACGT\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' \
                    > "bins/${sid}__$1_R1.fastq"
                if [ "$PAIRED" = "1" ]; then
                    printf '@r1\nACGTACGTACGTACGTACGTACGTACGT\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' \
                        > "bins/${sid}__$1_R2.fastq"
                fi
                echo "[STUB]   bin $1"
            }
            if ~{if extract_unclassified then "true" else "false"}; then emit_bin unclassified; fi
            for t in ~{sep=' ' extract_taxa_denovo}; do emit_bin "taxon_$t"; done
            for t in ~{sep=' ' extract_taxa_ref};    do emit_bin "ref_$t";   done
        done
'''

STUB_CONTIGS = r'''
        set -euo pipefail
        echo "[STUB] ClassifyContigs — database not localized"
        mkdir -p reports
        printf 'assembly\tclade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n' > reports/contig_reports.tsv
        for f in ~{sep=' ' assemblies}; do
            [ -s "$f" ] || continue
            b=$(basename "$f" .assembly.fasta)
            printf '%s\t100.0000\t1\t1\tspecies\t562\tEscherichia coli\n' "$b" >> reports/contig_reports.tsv
        done
'''

# swap the heavy containers for a tiny one — the stubs are pure shell
STUB_DOCKER = 'debian:stable-slim'


def replace_command(text, task_name, new_body, new_docker=None):
    """Replace the command <<< ... >>> body of a named task."""
    tstart = text.index('task %s {' % task_name)
    tend = text.index('\ntask ', tstart + 1) if '\ntask ' in text[tstart + 1:] else len(text)
    tend = text.index('\ntask ', tstart + 1) if text.find('\ntask ', tstart + 1) != -1 else len(text)
    block = text[tstart:tend]

    cstart = block.index('command <<<')
    cend = block.index('>>>', cstart)
    newblock = block[:cstart] + 'command <<<\n' + new_body + '\n    ' + block[cend:]

    if new_docker:
        newblock = re.sub(r'docker: "[^"]+"', 'docker: "%s"' % new_docker, newblock)
    return text[:tstart] + newblock + text[tend:]


def main():
    if not os.path.exists(WDL):
        sys.exit("run build_wdl.py first")

    # the stub never reads the database, but WDL still requires the File to
    # exist, so make the fixture self-contained
    work = os.path.join(HERE, 'work')
    os.makedirs(work, exist_ok=True)
    dummy = os.path.join(work, 'dummy_db.tar.gz')
    if not os.path.exists(dummy):
        import gzip
        with gzip.open(dummy, 'wb') as fh:
            fh.write(b'stub')
        print("created", dummy)

    with open(WDL) as fh:
        text = fh.read()

    text = replace_command(text, 'MetabuliClassifyAndExtract', STUB_CLASSIFY, STUB_DOCKER)
    text = replace_command(text, 'ClassifyContigs', STUB_CONTIGS, STUB_DOCKER)
    text = text.replace(
        '## ont-metabuli — Terra/WDL port',
        '## ont-metabuli — STUBBED COPY FOR LOCAL TESTING (do not run on Terra)')

    with open(OUT, 'w') as fh:
        fh.write(text)
    print("wrote", OUT)


if __name__ == '__main__':
    main()
