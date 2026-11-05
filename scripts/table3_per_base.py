#!/usr/bin/env python3
"""Manuscript Table 3: per-base fraction of reads called 8-oxoG at the designed positions.

Input: <results>/<construct>.filt.bed from run_plasmid_inference.sh
       (modkit pileup with --filter-threshold G:0.9980469).

Writes two TSVs:
  table3_per_base.tsv  one row per 8-oxoG construct at ITS OWN designed position + control value there
  table3_grid.tsv      every construct (incl. control) at every designed position (transparency)
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_samples, control_name, read_bedmethyl  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="output folder of run_plasmid_inference.sh")
    ap.add_argument("--samples", default=os.path.join(os.path.dirname(__file__), "..", "samples.tsv"))
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    samples = read_samples(a.samples)
    ctl = control_name(samples)
    oxo = [r for r in samples if r["designed_pos0"] is not None]
    data = {r["construct"]: read_bedmethyl(os.path.join(a.results, f"{r['construct']}.filt.bed"), plus_only=True)
            for r in samples}
    os.makedirs(a.outdir, exist_ok=True)

    def val(construct, pos):
        v = data[construct].get(pos)
        return (f"{v[0]:.1f}", str(v[1])) if v else ("NA", "0")

    out1 = os.path.join(a.outdir, "table3_per_base.tsv")
    with open(out1, "w") as fh:
        fh.write("construct\tinsert_offset\tpos0\tpct_8oxoG\tn_valid\tcontrol_pct_8oxoG\tcontrol_n_valid\n")
        for r in oxo:
            pct, n = val(r["construct"], r["designed_pos0"])
            cpct, cn = val(ctl, r["designed_pos0"])
            fh.write(f"{r['construct']}\t{r['insert_offset']}\t{r['designed_pos0']}\t{pct}\t{n}\t{cpct}\t{cn}\n")

    out2 = os.path.join(a.outdir, "table3_grid.tsv")
    with open(out2, "w") as fh:
        fh.write("pos0\tinsert_offset\t" + "\t".join(r["construct"] for r in samples) + "\n")
        for r in oxo:
            cells = [val(s["construct"], r["designed_pos0"])[0] for s in samples]
            fh.write(f"{r['designed_pos0']}\t{r['insert_offset']}\t" + "\t".join(cells) + "\n")

    for path in (out1, out2):
        print(f"--- {path}")
        print(open(path).read())


if __name__ == "__main__":
    main()
