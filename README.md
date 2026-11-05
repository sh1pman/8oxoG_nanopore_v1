# 8-oxoG nanopore decoder v1 (ONT R10.4.1)

Trained model and the minimal data needed to **reproduce the inference results** of:

> Sergeev A.V., Genatullina A.I., Grigorieva M.A., Chuprikova M.M., Meschaninova M.I., Zharkov D.O., Zvereva M.E.
> *DNA libraries for training a neural network model for nanopore signal decoding: the case of oxidative damage to guanosine residues* (translated title).
> Biomedical Chemistry: Research and Methods, 2026 (submitted).

The model is a Remora ConvLSTM (`ConvLSTM_w_ref`) binary classifier, canonical G vs 8-oxo-dG
(modification code `o`, ChEBI:44605), exported for Dorado. This repository does not cover training;
it lets you re-run the model on the validation data and obtain the numbers reported in the paper.

## What you can reproduce

| Result in the paper | Script | Data |
|---|---|---|
| Table 3: per-base 8-oxoG calls at the designed positions of the plasmid constructs | `scripts/table3_per_base.py` | plasmid POD5 (Zenodo) |
| Table 4 and signal decay: sliding-window analysis vs. control and plasmid background | `scripts/table4_window_scan.py` | plasmid POD5 (Zenodo) |
| 99.9% per-chunk accuracy on held-out cassette data (canonical G vs 8-oxoG, same site) | `scripts/validate_ext_val.sh` | `ext_val/` (this repo) |

## Contents

```
model/
  dorado_8oxoG_model_v1/      exported model for dorado --modified-bases-models
  remora_8oxoG_v1.pt          Remora checkpoint (used for the published numbers)
  SHA256SUMS
reference/
  plasmid.fasta               record "pUC28-rc": pUC19 (2686 bp), reverse-complement orientation
  insert.fasta                24-nt insert region, found at 0-based 2585-2608 of plasmid.fasta
  designed_sites.bed          the four 8-oxoG positions (0-based, + strand)
samples.tsv                   construct <-> barcode <-> designed position
ext_val/
  ext_val_neg/                Remora chunks, lib_G, target site, canonical G (2,217 chunks)
  ext_val_pos/                Remora chunks, lib_oG, target site, 8-oxoG     (2,044 chunks)
demo/pod5/                    500 reads per construct (random, seed 42) for a quick test
expected/                     outputs of the full and demo runs, for comparison
scripts/                      pipeline and table scripts
```

## Constructs

Linearised plasmids with a single 8-oxoG introduced by ligation of a chemically synthesised,
modified oligonucleotide into an nCas9-generated gap (no cloning or passage through cells, so
the lesion is expected at ~100% of molecules). Sequenced with SQK-NBD114-24 on R10.4.1 (5 kHz).

| construct | barcode | 8-oxoG position, 0-based in `pUC28-rc` | position in insert |
|---|---|---|---|
| 38G (control, no 8-oxoG) | barcode17 | - | - |
| 38oxoG1 | barcode18 | 2587 | 2 |
| 38oxoG2 | barcode19 | 2592 | 7 |
| 38oxoG3 | barcode21 | 2599 | 14 |
| 38oxoG4 | barcode20 | 2604 | 19 |

Note barcode21 = 38oxoG3 and barcode20 = 38oxoG4 (not in numerical order).

## Requirements

| Software | Version used |
|---|---|
| Dorado | 1.1.1 (>= 1.0 required), basecaller model `sup@v5.2.0` = `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` |
| Remora (ont-remora) | 3.3.0 |
| modkit | 0.5.0 |
| samtools | 1.17-1.18 (any recent version) |
| Python | 3.10 (standard library only for the table scripts) |

The original plasmid basecalling used Dorado 1.1.1 on an NVIDIA RTX 3090; Remora inference used an
RTX 3060 Ti (8 GB), Ubuntu under WSL2. Remora also runs on CPU (`DEVICE=cpu`), slowly.

Dorado options for the scripts: `DORADO=/path/to/dorado` if 1.1.1 is not the `dorado` in your PATH;
`MODELS_DIR=...` for the basecaller model folder (default `~/.cache/dorado_models`, downloaded on
first use).

```bash
conda env create -f environment.yml && conda activate oxog-v1
# Dorado 1.1.1 and modkit 0.5.0 are standalone binaries:
#   https://github.com/nanoporetech/dorado   https://github.com/nanoporetech/modkit/releases/tag/v0.5.0
```

## Quick test with the demo reads (minutes)

```bash
git clone https://github.com/sh1pman/8oxoG_nanopore_v1 && cd 8oxoG_nanopore_v1
POD5_SUFFIX=_demo.pod5 bash scripts/reproduce_all.sh demo/pod5 results_demo
```

500 reads per construct check that the whole pipeline runs on your system; compare with
`expected/demo/`. With this few reads only the two strongest constructs stand out in the window
analysis (38oxoG2 and 38oxoG3: 100th percentile, 3-5x the background p99), while 38oxoG1 and
38oxoG4 stay within background noise (91-99.6th percentile). The published values need the full data.

## Full reproduction

Download the plasmid reads (one POD5 per construct, ~2 GB) from Zenodo: **DOI: 10.5281/zenodo.23977043**, then

```bash
sha256sum -c SHA256SUMS                       # inside the downloaded folder
bash scripts/reproduce_all.sh /path/to/pod5 results
```

Outputs: `results/tables/*.tsv` (compare with `expected/full/`), `results/ext_val/*.txt`.

### Expected results

Table 3 - per-base fraction of reads called 8-oxoG at the designed position (`table3_per_base.tsv`):

| construct | % 8-oxoG at its designed position | control 38G at the same position, % |
|---|---|---|
| 38oxoG1 | 1.7 | 0.2 |
| 38oxoG2 | 5.4 | 0.1 |
| 38oxoG3 | 26.7 | 2.1 |
| 38oxoG4 | 0.0 | 0.0 |

Table 4 - window centred on each construct's own designed position, construct minus control
(`table4_own_position.tsv`; background = all windows outside the insert):

| construct | delta, +/-12 nt | percentile | delta / background p99 | delta, +/-5 nt | percentile |
|---|---|---|---|---|---|
| 38oxoG1 | 1.47 | 100 | 4.9 | 2.20 | 100 |
| 38oxoG2 | 3.54 | 100 | 10.2 | 5.82 | 100 |
| 38oxoG3 | 5.24 | 100 | 17.2 | 9.06 | 100 |
| 38oxoG4 | 1.89 | 100 | 6.1 | 2.30 | 100 |

Cassette held-out sets (canonical G vs 8-oxoG at the same site and anchor context):

| file | chunks | correct | accuracy | after filtering the 10% least confident calls |
|---|---|---|---|---|
| `ext_val_neg_all.txt` (canonical G) | 2,217 | 2,211 | 99.73% (specificity) | 100% |
| `ext_val_pos_all.txt` (8-oxoG) | 2,044 | 2,039 | 99.76% (sensitivity) | 100% |
| both | 4,261 | 4,250 | 99.74% | 100% |
| `ext_val_balanced.txt` (balanced batch) | 2,048 | 2,046 | 99.90% | 100% |

The model checkpoint was selected on the balanced 2,048-chunk batch during training, so 99.9% is
slightly optimistic; the all-chunk figures are the conservative estimate. `remora validate` splits
the balanced batch between the two sets at random, so its confusion matrix varies between runs
(e.g. [[1022, 0], [2, 1024]]); the per-set results are deterministic.

These values come from re-running this pipeline from the POD5 files (Dorado 1.1.1 on an
RTX 3060 Ti). Basecalling on a different GPU model can let slightly different reads pass alignment
and barcode filtering; with the original basecalls (RTX 3090) the values differed by at most 0.1
(e.g. 26.8 vs 26.7 in Table 3). Conclusions are unaffected.

## How the numbers are computed (please read before comparing)

1. **Inference is reference-anchored Remora, not Dorado.** `dorado --modified-bases-models` scores
   only positions basecalled as G. 8-oxoG is often basecalled as another base, so Dorado reports
   lower fractions at lesion sites. All published numbers use
   `remora infer from_pod5_and_bam --reference-anchored` with `model/remora_8oxoG_v1.pt`.
2. **Table 3 uses a fixed filter threshold.** modkit normally estimates a confidence threshold
   per run; the original run gave 0.9980469 for every barcode, so it is pinned
   (`--filter-threshold G:0.9980469`) to make the result independent of the run.
3. **Table 4 uses unfiltered calls on the modified strand only.** `modkit pileup --no-filtering`;
   only `+`-strand (reference-G) rows, because 8-oxoG is on that strand. A window centred at c
   covers positions c-12 ... c+12 (25 positions; "24-nt window" in the paper) or c-5 ... c+5.
   Its value is the mean % modified over the `+`-strand G positions inside it.
4. **Compare each construct at its own position.** The four designed positions lie within 17 nt
   of each other, so a +/-12 window centred on one designed position also contains the others;
   a strong lesion raises the neighbouring windows too. `table4_grid.tsv` shows every
   construct x every window centre to make this visible.
5. **Per-base calls are spread to neighbouring bases.** The model detects the lesion within a few
   nucleotides but often not on the exact base (e.g. 38oxoG4: 0% at its own base, yet its window is
   above every background window). Signal decays to background within ~4-6 nt of the insert
   (`decay_<construct>.tsv`).

## Using the model on new data

```bash
# convenient, basecall-anchored (see note 1 above)
dorado basecaller sup@v5.2.0 reads.pod5 --emit-moves \
  --modified-bases-models model/dorado_8oxoG_model_v1 --reference ref.fa > calls.bam

# reference-anchored, as in the paper
dorado basecaller sup@v5.2.0 reads.pod5 --emit-moves --reference ref.fa > mapped.bam
samtools sort -o mapped.sorted.bam mapped.bam && samtools index mapped.sorted.bam
remora infer from_pod5_and_bam reads.pod5 mapped.sorted.bam --model model/remora_8oxoG_v1.pt \
  --reference-anchored --out-bam calls.bam --device 0
```

Model settings: motif `G` offset 0; chunk context 50/50 samples; k-mer context 4/4 (9-mer);
rough signal rescaling with the R10.4.1 9-mer level table (embedded in both model files).

Limitations: trained on synthetic DNA cassettes plus genomic amplicon negatives; localisation is
~10 nt rather than single-base; some G-run contexts give elevated background calls on canonical G.
Use matched unmodified controls when interpreting calls.

## License and citation

Code and model: MIT (see `LICENSE`). Data on Zenodo: CC BY 4.0.
Please cite the paper above; machine-readable metadata in `CITATION.cff`.
Contact: A.V. Sergeev, sergeevav_01@my.msu.ru
