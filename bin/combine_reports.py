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
