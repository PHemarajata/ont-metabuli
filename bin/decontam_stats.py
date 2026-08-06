#!/usr/bin/env python3
"""
Significance / decontamination statistics for ONT metagenomics.

Two modes, auto-selected from the samplesheet:

  (A) NEGATIVE-CONTROL mode  — when >=1 sample is flagged as a control.
      For each biological sample and each taxon, compare the taxon's proportion
      of classified reads in the sample vs. in the (batch-matched, else pooled)
      negative control with a one-sided Fisher exact test (enriched in sample).
      A taxon is called:
        contaminant       if its proportion in the control >= in the sample
        significant       if BH q < alpha AND fold-enrichment >= min_fold
        below_threshold   if reads < min_reads or rel-abundance < min_abundance
        not_significant   otherwise
      This is the frequency-based logic behind decontam-style contaminant
      identification, adapted for a single control per batch.

  (B) NO-CONTROL mode  — when no control is present (the common ONT case).
      Significance testing without a blank is inherently limited: there is no
      empirical contamination background to test against. We therefore report a
      transparent HEURISTIC: for each taxon, a one-sided Poisson tail test of the
      observed read count against an assumed low background rate (--bg-rate x
      classified reads), plus a Wilson 95% CI on relative abundance. A taxon is:
        detected          if reads>=min_reads, rel-abundance>=min_abundance, BH q<alpha
        below_threshold   if it fails the count/abundance floor
        not_significant   otherwise
      Treat these as detection confidence flags, NOT contamination-corrected
      calls. Include a negative control whenever possible for rigorous results.

Outputs: taxon_stats.tsv, significant_taxa.tsv, sankey_flags.tsv, stats_summary.json
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, poisson

CONTROL_CALL = 'contaminant'


def bh_qvalues(pvals):
    """Benjamini-Hochberg FDR. Returns array aligned to input order."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def wilson_ci(count, total, z=1.96):
    if total == 0:
        return (0.0, 0.0)
    phat = count / total
    denom = 1 + z**2 / total
    centre = (phat + z**2 / (2 * total)) / denom
    half = (z * math.sqrt(phat * (1 - phat) / total + z**2 / (4 * total**2))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--abundance', required=True)   # combined_abundance_long.tsv
    ap.add_argument('--metadata', required=True)    # sample sample_type batch
    ap.add_argument('--rank', default='species')
    ap.add_argument('--min-reads', type=int, default=10)
    ap.add_argument('--min-abundance', type=float, default=1e-4)
    ap.add_argument('--alpha', type=float, default=0.05)
    ap.add_argument('--min-fold', type=float, default=5.0)
    ap.add_argument('--bg-rate', type=float, default=1e-5)
    ap.add_argument('--out-prefix', default='.')
    args = ap.parse_args()
    os.makedirs(args.out_prefix, exist_ok=True)
    eps = 1e-9

    abund = pd.read_csv(args.abundance, sep='\t', dtype={'taxid': str})
    meta = pd.read_csv(args.metadata, sep='\t', dtype=str).fillna('all')
    meta['sample_type'] = meta['sample_type'].str.lower()
    controls = meta.loc[meta['sample_type'] == 'control', 'sample'].tolist()
    bio_samples = meta.loc[meta['sample_type'] != 'control', 'sample'].tolist()
    batch_of = dict(zip(meta['sample'], meta['batch']))

    # per-sample classified total (constant within a sample)
    totals = abund.groupby('sample')['total_classified'].max().to_dict()
    # reads[sample][taxid] = clade_reads ; names[taxid] = name
    reads = {}
    names = {}
    for _, r in abund.iterrows():
        reads.setdefault(r['sample'], {})[r['taxid']] = int(r['clade_reads'])
        names[r['taxid']] = r['name']

    mode = 'negative_control' if controls else 'no_control'
    out_rows = []
    per_sample_sig = {}

    def control_pool(sample):
        """Return (reads_dict, total) for controls matched to this sample."""
        b = batch_of.get(sample, 'all')
        matched = [c for c in controls if batch_of.get(c, 'all') == b]
        if not matched:
            matched = controls  # fall back to all controls pooled
        pooled = {}
        tot = 0
        for c in matched:
            tot += totals.get(c, 0)
            for tx, ct in reads.get(c, {}).items():
                pooled[tx] = pooled.get(tx, 0) + ct
        return pooled, tot

    for sample in bio_samples:
        A = totals.get(sample, 0)
        if A == 0:
            continue
        s_reads = reads.get(sample, {})
        taxa = list(s_reads.keys())

        if mode == 'negative_control':
            c_reads, C = control_pool(sample)
            pvals, recs = [], []
            for tx in taxa:
                a = s_reads[tx]
                c = c_reads.get(tx, 0)
                prop_s = a / A
                prop_c = (c / C) if C else 0.0
                table = [[a, max(A - a, 0)], [c, max(C - c, 0)]]
                try:
                    _, p = fisher_exact(table, alternative='greater')
                except ValueError:
                    p = 1.0
                fold = (prop_s + eps) / (prop_c + eps)
                pvals.append(p)
                recs.append((tx, a, prop_s, c, prop_c, fold, p))
            qvals = bh_qvalues(pvals)
            nsig = 0
            for (tx, a, prop_s, c, prop_c, fold, p), q in zip(recs, qvals):
                if a < args.min_reads or prop_s < args.min_abundance:
                    call = 'below_threshold'
                elif prop_c >= prop_s:
                    call = CONTROL_CALL
                elif q < args.alpha and fold >= args.min_fold:
                    call = 'significant'
                else:
                    call = 'not_significant'
                if call == 'significant':
                    nsig += 1
                lo, hi = wilson_ci(a, A)
                out_rows.append([sample, tx, names.get(tx, tx), args.rank, a, A,
                                 f'{prop_s:.6g}', f'{lo:.4g}', f'{hi:.4g}',
                                 c, f'{prop_c:.6g}', f'{fold:.4g}',
                                 f'{p:.4g}', f'{q:.4g}', call, mode])
            per_sample_sig[sample] = nsig

        else:  # no_control
            pvals, recs = [], []
            for tx in taxa:
                a = s_reads[tx]
                prop_s = a / A
                lam = max(args.bg_rate * A, eps)
                p = float(poisson.sf(a - 1, lam))  # P(X >= a)
                pvals.append(p)
                recs.append((tx, a, prop_s, p))
            qvals = bh_qvalues(pvals)
            nsig = 0
            for (tx, a, prop_s, p), q in zip(recs, qvals):
                if a < args.min_reads or prop_s < args.min_abundance:
                    call = 'below_threshold'
                elif q < args.alpha:
                    call = 'detected'
                else:
                    call = 'not_significant'
                if call == 'detected':
                    nsig += 1
                lo, hi = wilson_ci(a, A)
                out_rows.append([sample, tx, names.get(tx, tx), args.rank, a, A,
                                 f'{prop_s:.6g}', f'{lo:.4g}', f'{hi:.4g}',
                                 '', '', '', f'{p:.4g}', f'{q:.4g}', call, mode])
            per_sample_sig[sample] = nsig

    cols = ['sample', 'taxid', 'name', 'rank', 'reads', 'total_classified',
            'rel_abundance', 'ci_low', 'ci_high', 'control_reads',
            'control_rel_abundance', 'fold_enrichment', 'p_value', 'q_value',
            'call', 'mode']
    df = pd.DataFrame(out_rows, columns=cols)
    df = df.sort_values(['sample', 'reads'], ascending=[True, False])

    p = lambda f: os.path.join(args.out_prefix, f)
    df.to_csv(p('taxon_stats.tsv'), sep='\t', index=False)

    sig_calls = {'significant', 'detected'}
    df[df['call'].isin(sig_calls)].to_csv(p('significant_taxa.tsv'), sep='\t', index=False)

    df[['sample', 'taxid', 'call']].to_csv(p('sankey_flags.tsv'), sep='\t', index=False)

    summary = {
        'mode': mode,
        'n_biological_samples': len(bio_samples),
        'n_controls': len(controls),
        'controls': controls,
        'rank': args.rank,
        'params': {'min_reads': args.min_reads, 'min_abundance': args.min_abundance,
                   'alpha': args.alpha, 'min_fold': args.min_fold,
                   'bg_rate': args.bg_rate},
        'significant_per_sample': per_sample_sig,
        'n_significant_total': int(sum(per_sample_sig.values())),
        'caveat': ('No negative control supplied: calls are Poisson-based '
                   'detection-confidence heuristics, not contamination-corrected.'
                   if mode == 'no_control' else
                   'Contaminant calls are frequency-based vs. batch-matched controls.')
    }
    with open(p('stats_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2)

    print(f"[decontam_stats] mode={mode} controls={len(controls)} "
          f"sig_total={summary['n_significant_total']}", file=sys.stderr)


if __name__ == '__main__':
    main()
