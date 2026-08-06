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
