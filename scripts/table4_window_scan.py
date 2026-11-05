#!/usr/bin/env python3
"""Manuscript Table 4 (+ signal decay): sliding-window analysis on the plasmid constructs.

Input: <results>/<construct>.nofilt.bed from run_plasmid_inference.sh (modkit pileup --no-filtering).

Definitions (identical to the analysis in the manuscript):
  * Only '+' strand rows are used, i.e. reference-G positions - the strand that carries 8-oxoG.
  * A window centred at c with half-width h spans reference positions c-h .. c+h inclusive
    (h = 12 -> 25 positions, reported as the "24 nt" window; h = 5 -> 11 positions, "+/-5 nt").
  * Window value = mean percent_modified over the '+'-strand G positions inside the window.
  * delta = window value (construct) - window value (control 38G), same window.
  * Background null = deltas of every window lying entirely outside the insert.
    Percentile = % of background windows with delta <= the window of interest.

Writes:
  table4_own_position.tsv  each 8-oxoG construct, window centred on ITS OWN designed position
  table4_grid.tsv          each 8-oxoG construct x each designed position (shows window overlap)
  table4_background.tsv    background-null summary per construct and window size
  decay_<construct>.tsv    window walked across the insert +/-48 nt (step 2) for every construct
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_samples, control_name, read_fasta_seq, read_bedmethyl  # noqa: E402

HALF_WIDTHS = (12, 5)      # 24-nt window (+/-12) and +/-5 nt window
DECAY_SPAN, DECAY_STEP = 48, 2


def wmean(d, lo, hi):
    vals = [d[p][0] for p in range(lo, hi + 1) if p in d]
    return sum(vals) / len(vals) if vals else None


def pctile(v, dist):
    return 100.0 * sum(x <= v for x in dist) / len(dist)


def q(sorted_vals, frac):
    n = len(sorted_vals)
    return sorted_vals[min(n - 1, int(frac * n))]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--samples", default=os.path.join(here, "..", "samples.tsv"))
    ap.add_argument("--reference", default=os.path.join(here, "..", "reference", "plasmid.fasta"))
    ap.add_argument("--insert", default=os.path.join(here, "..", "reference", "insert.fasta"))
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    _, plasmid = read_fasta_seq(a.reference)
    _, insert = read_fasta_seq(a.insert)
    ins_start = plasmid.find(insert)
    if ins_start < 0:
        raise SystemExit("insert sequence not found (forward strand) in the plasmid reference")
    ins_end = ins_start + len(insert) - 1
    plen = len(plasmid)
    print(f"plasmid length {plen}; insert at 0-based {ins_start}-{ins_end}")

    samples = read_samples(a.samples)
    ctl = control_name(samples)
    oxo = [r for r in samples if r["designed_pos0"] is not None]
    data = {r["construct"]: read_bedmethyl(os.path.join(a.results, f"{r['construct']}.nofilt.bed"), plus_only=True)
            for r in samples}
    os.makedirs(a.outdir, exist_ok=True)

    def delta_at(construct, center, half):
        lo, hi = center - half, center + half
        cm, om = wmean(data[ctl], lo, hi), wmean(data[construct], lo, hi)
        if cm is None or om is None:
            return None, om, cm
        return om - cm, om, cm

    # background null per construct and half-width
    null = {}
    for half in HALF_WIDTHS:
        for r in oxo:
            vals = []
            for center in range(half, plen - half):
                lo, hi = center - half, center + half
                if hi < ins_start or lo > ins_end:
                    d, _, _ = delta_at(r["construct"], center, half)
                    if d is not None:
                        vals.append(d)
            null[(r["construct"], half)] = sorted(vals)

    with open(os.path.join(a.outdir, "table4_background.tsv"), "w") as fh:
        fh.write("construct\thalf_width\tn_windows\tmean\tp95\tp99\tmax\n")
        for (c, half), dd in null.items():
            fh.write(f"{c}\t{half}\t{len(dd)}\t{sum(dd)/len(dd):.3f}\t{q(dd,.95):.2f}\t{q(dd,.99):.2f}\t{max(dd):.2f}\n")

    with open(os.path.join(a.outdir, "table4_own_position.tsv"), "w") as fh:
        fh.write("construct\tinsert_offset\tpos0\thalf_width\tox_window_pct\tctl_window_pct\tdelta\t"
                 "percentile_vs_background\tbackground_p99\tdelta_over_p99\n")
        for r in oxo:
            for half in HALF_WIDTHS:
                d, om, cm = delta_at(r["construct"], r["designed_pos0"], half)
                dd = null[(r["construct"], half)]
                p99 = q(dd, .99)
                fh.write(f"{r['construct']}\t{r['insert_offset']}\t{r['designed_pos0']}\t{half}\t{om:.2f}\t{cm:.2f}\t"
                         f"{d:.2f}\t{pctile(d, dd):.1f}\t{p99:.2f}\t{d/p99:.1f}\n")

    with open(os.path.join(a.outdir, "table4_grid.tsv"), "w") as fh:
        fh.write("half_width\twindow_centre_pos0\twindow_centre_construct\t" +
                 "\t".join(f"delta_{r['construct']}" for r in oxo) + "\n")
        for half in HALF_WIDTHS:
            for site in oxo:
                cells = []
                for r in oxo:
                    d, _, _ = delta_at(r["construct"], site["designed_pos0"], half)
                    cells.append(f"{d:.2f}" if d is not None else "NA")
                fh.write(f"{half}\t{site['designed_pos0']}\t{site['construct']}\t" + "\t".join(cells) + "\n")

    half = 12
    ins_center = (ins_start + ins_end) // 2
    for r in oxo:
        dd = null[(r["construct"], half)]
        with open(os.path.join(a.outdir, f"decay_{r['construct']}.tsv"), "w") as fh:
            fh.write("window_centre_pos0\tdist_from_insert_edge\twindow_overlaps_insert\tdelta\tpercentile_vs_background\n")
            for center in range(ins_center - DECAY_SPAN - half, ins_center + DECAY_SPAN + half + 1, DECAY_STEP):
                d, _, _ = delta_at(r["construct"], center, half)
                if d is None:
                    continue
                dist = ins_start - center if center < ins_start else (center - ins_end if center > ins_end else 0)
                overlap = (ins_start - half) <= center <= (ins_end + half)
                fh.write(f"{center}\t{dist}\t{'yes' if overlap else 'no'}\t{d:.2f}\t{pctile(d, dd):.1f}\n")

    for name in ("table4_own_position.tsv", "table4_background.tsv"):
        path = os.path.join(a.outdir, name)
        print(f"--- {path}")
        print(open(path).read())


if __name__ == "__main__":
    main()
