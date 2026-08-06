#!/usr/bin/env python3
"""
Build a self-contained interactive taxonomic Sankey (Domain -> ... -> Species)
from one Metabuli/Kraken2-style report.tsv, using plotly.

Flow width = clade reads. Nodes are laid out in rank columns. Species nodes are
recoloured by their decontam call (from sankey_flags.tsv) when available:
significant/detected = green, contaminant = grey, below/ns = muted.

The output HTML embeds plotly.js (include_plotlyjs=True) so it works offline.
"""
import argparse
import os
import sys

# canonical single-letter rank codes we may place as columns
RANK_NAMES = {'D': 'Domain', 'K': 'Kingdom', 'P': 'Phylum', 'C': 'Class',
              'O': 'Order', 'F': 'Family', 'G': 'Genus', 'S': 'Species'}

# Metabuli v1.2.0 uses full rank names; older format uses single letters.
_FULL2LETTER = {'superkingdom': 'D', 'domain': 'D', 'kingdom': 'K', 'phylum': 'P',
                'class': 'C', 'order': 'O', 'family': 'F', 'genus': 'G', 'species': 'S'}


def norm_rank(rank):
    r = (rank or '').strip()
    if len(r) == 1 and r.upper() in 'DKPCOFGS':
        return r.upper()
    return _FULL2LETTER.get(r.lower())

DOMAIN_COLORS = {
    'Bacteria': '#4C78A8', 'Viruses': '#E45756', 'Archaea': '#B279A2',
    'Eukaryota': '#F58518', 'Fungi': '#72B7B2', 'unknown': '#9D9D9D',
}
CALL_COLORS = {
    'significant': '#2CA02C', 'detected': '#2CA02C',
    'contaminant': '#8C8C8C', 'below_threshold': '#D9D9D9',
    'not_significant': '#BFBFBF',
}


def parse_tree(path):
    """Parse report into ordered nodes with parent pointers (via depth stack)."""
    nodes = []
    stack = {}  # depth -> node index
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            try:
                clade_reads = int(parts[1])
            except ValueError:
                continue
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name_raw = parts[5]
            depth = (len(name_raw) - len(name_raw.lstrip(' '))) // 2
            name = name_raw.strip()
            if taxid == '0':  # unclassified
                continue
            parent = None
            for d in range(depth - 1, -1, -1):
                if d in stack:
                    parent = stack[d]
                    break
            idx = len(nodes)
            nodes.append({'taxid': taxid, 'rank': rank, 'name': name,
                          'clade_reads': clade_reads, 'depth': depth,
                          'parent': parent})
            stack[depth] = idx
            # drop deeper stack entries
            for d in list(stack):
                if d > depth:
                    del stack[d]
    return nodes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--report', required=True)
    ap.add_argument('--sample', required=True)
    ap.add_argument('--flags', default=None)
    ap.add_argument('--ranks', default='D,P,C,O,F,G,S')
    ap.add_argument('--top-n', type=int, default=15)
    ap.add_argument('--min-reads', type=int, default=5)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    rank_order = [r.strip().upper() for r in args.ranks.split(',') if r.strip()]
    rank_pos = {r: i for i, r in enumerate(rank_order)}

    # ---- flags for this sample ----
    call_of = {}
    if args.flags and os.path.exists(args.flags):
        with open(args.flags) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            try:
                si, ti, ci = header.index('sample'), header.index('taxid'), header.index('call')
            except ValueError:
                si = ti = ci = None
            if si is not None:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    if len(f) > max(si, ti, ci) and f[si] == args.sample:
                        call_of[f[ti]] = f[ci]

    nodes = parse_tree(args.report)

    def main_rank(code):
        l = norm_rank(code)
        return l if (l and l in rank_pos) else None

    # canonical candidate nodes (rank in requested set, >= min reads)
    canon = [i for i, n in enumerate(nodes)
             if main_rank(n['rank']) and n['clade_reads'] >= args.min_reads]

    # keep top-N per rank
    keep = set()
    by_rank = {}
    for i in canon:
        by_rank.setdefault(main_rank(nodes[i]['rank']), []).append(i)
    for r, idxs in by_rank.items():
        idxs.sort(key=lambda i: -nodes[i]['clade_reads'])
        keep.update(idxs[:args.top_n])

    def nearest_kept_ancestor(i):
        p = nodes[i]['parent']
        while p is not None:
            if p in keep:
                return p
            p = nodes[p]['parent']
        return None

    def domain_of(i):
        p = i
        seen = 0
        while p is not None and seen < 50:
            if norm_rank(nodes[p]['rank']) == 'D':
                return nodes[p]['name']
            p = nodes[p]['parent']
            seen += 1
        return 'unknown'

    # collapse to only the ranks actually present, so columns fill full width
    present = sorted({main_rank(nodes[i]['rank']) for i in keep}, key=lambda r: rank_pos[r])
    col_of = {r: k for k, r in enumerate(present)}
    ncol = max(len(present) - 1, 1)

    # ---- build plotly node/link arrays ----
    node_ids = {}         # keep-index -> plotly node index
    labels, node_colors, node_x = [], [], []
    for i in sorted(keep, key=lambda i: (rank_pos[main_rank(nodes[i]['rank'])], -nodes[i]['clade_reads'])):
        node_ids[i] = len(labels)
        n = nodes[i]
        rm = main_rank(n['rank'])
        labels.append(f"{n['name']} ({n['clade_reads']:,})")
        # colour: species by call if available, else by domain
        if rm == 'S' and n['taxid'] in call_of:
            node_colors.append(CALL_COLORS.get(call_of[n['taxid']], '#BFBFBF'))
        else:
            node_colors.append(DOMAIN_COLORS.get(domain_of(i), DOMAIN_COLORS['unknown']))
        node_x.append(0.01 + 0.98 * col_of[rm] / ncol)

    src, tgt, val, lcolor = [], [], [], []
    for i in keep:
        anc = nearest_kept_ancestor(i)
        if anc is None:
            continue
        src.append(node_ids[anc])
        tgt.append(node_ids[i])
        val.append(max(nodes[i]['clade_reads'], 1))
        base = DOMAIN_COLORS.get(domain_of(i), DOMAIN_COLORS['unknown'])
        lcolor.append(_rgba(base, 0.35))

    _render(args, labels, node_colors, node_x, src, tgt, val, lcolor, len(keep))
    print(f"[make_sankey] {args.sample}: {len(keep)} nodes, {len(src)} links",
          file=sys.stderr)


def _rgba(hexc, a):
    hexc = hexc.lstrip('#')
    r, g, b = int(hexc[0:2], 16), int(hexc[2:4], 16), int(hexc[4:6], 16)
    return f"rgba({r},{g},{b},{a})"


def _render(args, labels, node_colors, node_x, src, tgt, val, lcolor, nkeep):
    title = f"Taxonomic flow — {args.sample}"
    if nkeep == 0:
        with open(args.out, 'w') as fh:
            fh.write(f"<html><body style='font-family:sans-serif'>"
                     f"<h3>{title}</h3><p>No taxa above the reporting "
                     f"threshold (min {args.min_reads} clade reads).</p>"
                     f"</body></html>")
        return
    import plotly.graph_objects as go
    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node=dict(label=labels, color=node_colors, x=node_x,
                  pad=12, thickness=16,
                  line=dict(color='rgba(0,0,0,0.25)', width=0.5),
                  hovertemplate='%{label}<extra></extra>'),
        link=dict(source=src, target=tgt, value=val, color=lcolor,
                  hovertemplate='%{source.label} → %{target.label}<br>'
                                '%{value:,} reads<extra></extra>'),
    ))
    fig.update_layout(
        title=dict(text=title, x=0.02, font=dict(size=18)),
        font=dict(size=11), margin=dict(l=10, r=10, t=50, b=30), height=720,
        annotations=[dict(text=("Green = statistically supported · Grey = "
                                "likely contaminant · Flow width = reads"),
                          showarrow=False, x=0.02, y=-0.06, xref='paper',
                          yref='paper', font=dict(size=10, color='#666'))])
    fig.write_html(args.out, include_plotlyjs=True, full_html=True)


if __name__ == '__main__':
    main()
