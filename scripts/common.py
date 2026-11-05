"""Shared helpers: read samples.tsv, FASTA and modkit bedMethyl files."""
import csv
import os


def read_samples(path):
    """Return list of dicts from samples.tsv (tab-separated, header row)."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["designed_pos0"] = None if r["designed_pos0"] == "NA" else int(r["designed_pos0"])
        r["insert_offset"] = None if r["insert_offset"] == "NA" else int(r["insert_offset"])
    return rows


def control_name(samples):
    ctl = [r["construct"] for r in samples if r["designed_pos0"] is None]
    if len(ctl) != 1:
        raise SystemExit("samples.tsv must contain exactly one control row (designed_pos0 = NA)")
    return ctl[0]


def read_fasta_seq(path):
    """Single-record FASTA -> (name, sequence)."""
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                seq.append(line.upper())
    return name, "".join(seq)


def read_bedmethyl(path, plus_only=False, mod_code="o"):
    """modkit bedMethyl -> {pos0: (percent_modified, n_valid_cov)}.

    Columns (modkit >= 0.3): 1 chrom, 2 start(0-based), 3 end, 4 mod code, 5 score, 6 strand,
    7-9 display, 10 N_valid_cov, 11 percent_modified, 12.. counts.
    Split on any whitespace so both tab- and space-delimited variants parse.
    """
    if not os.path.exists(path):
        raise SystemExit(f"missing pileup: {path}")
    d = {}
    with open(path) as fh:
        for line in fh:
            c = line.split()
            if len(c) < 11 or c[3] != mod_code:
                continue
            if plus_only and c[5] != "+":
                continue
            d[int(c[1])] = (float(c[10]), int(c[9]))
    return d
