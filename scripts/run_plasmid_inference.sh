#!/usr/bin/env bash
# Reproduce 8-oxoG calls on the plasmid validation constructs (manuscript Tables 3-4).
#
# For every construct in samples.tsv:
#   POD5 -> dorado basecaller (moves + inline alignment) -> samtools filter/sort
#        -> remora infer (reference-anchored, model v1) -> modkit pileup (unfiltered + fixed threshold)
#
# Usage:  bash scripts/run_plasmid_inference.sh <pod5_dir> <out_dir>
#   <pod5_dir>  folder with 38G.pod5, 38oxoG1.pod5 ... (Zenodo record), or demo/pod5 (files *_demo.pod5)
#   <out_dir>   created if missing
#
# Environment overrides:
#   DORADO          dorado executable (default: 'dorado' from PATH); must be >= 1.0, 1.1.1 was used
#   MODELS_DIR      where dorado looks for / downloads the basecaller model (default ~/.cache/dorado_models)
#   DEVICE=0        GPU id for Remora (DEVICE=cpu runs Remora on CPU; much slower)
#   BASECALL_MODEL  canonical basecaller model (default sup@v5.2.0 = dna_r10.4.1_e8.2_400bps_sup@v5.2.0)
#   BC_FILTER=1     keep only reads re-classified to the construct's barcode (0 = keep all reads)
#   POD5_SUFFIX     file suffix to look for (default ".pod5"; use "_demo.pod5" for the demo set)

set -euo pipefail

POD5_DIR=${1:?usage: run_plasmid_inference.sh <pod5_dir> <out_dir>}
OUT=${2:?usage: run_plasmid_inference.sh <pod5_dir> <out_dir>}

REPO=$(cd "$(dirname "$0")/.." && pwd)
REF="$REPO/reference/plasmid.fasta"
MODEL="$REPO/model/remora_8oxoG_v1.pt"
SAMPLES="$REPO/samples.tsv"

DORADO=${DORADO:-dorado}
MODELS_DIR=${MODELS_DIR:-$HOME/.cache/dorado_models}
# Model complex: dorado resolves it from the POD5 chemistry to dna_r10.4.1_e8.2_400bps_sup@v5.2.0
BASECALL_MODEL=${BASECALL_MODEL:-sup@v5.2.0}
KIT=SQK-NBD114-24
DEVICE=${DEVICE:-0}
BC_FILTER=${BC_FILTER:-1}
POD5_SUFFIX=${POD5_SUFFIX:-.pod5}
# Threshold that modkit estimated on the original run (identical for all barcodes).
# Pinned here so that the filtered pileup (Table 3) does not depend on per-run estimation.
FILTER_THRESHOLD=0.9980469

for tool in "$DORADO" samtools remora modkit; do
  command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found (check PATH or DORADO=...)" >&2; exit 1; }
done
# Run dorado with its own bundled libraries first, so a LD_LIBRARY_PATH left by another dorado
# installation (e.g. 0.9.x) cannot inject an incompatible libtorch.
DORADO_LIB="$(cd "$(dirname "$(command -v "$DORADO")")/.." && pwd)/lib"
run_dorado() {
  if [[ -d "$DORADO_LIB" ]]; then
    LD_LIBRARY_PATH="$DORADO_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$DORADO" "$@"
  else
    "$DORADO" "$@"
  fi
}
dorado_ver=$(run_dorado --version 2>&1 | tail -n 1 | tr -d '\r')
if [[ ! "$dorado_ver" =~ ^[1-9][0-9]*\. ]]; then
  echo "ERROR: $DORADO reports version '$dorado_ver'; dorado >= 1.0 is required (1.1.1 was used)." >&2
  echo "       Install 1.1.1 alongside and run with DORADO=/path/to/dorado-1.1.1-linux-x64/bin/dorado" >&2
  exit 1
fi
[[ "$dorado_ver" == 1.1.1* ]] || echo "WARNING: dorado $dorado_ver; published numbers were produced with 1.1.1" >&2
mkdir -p "$MODELS_DIR"
[[ -f "$MODEL" ]] || { echo "ERROR: model not found: $MODEL" >&2; exit 1; }
[[ -f "$REF" ]]   || { echo "ERROR: reference not found: $REF" >&2; exit 1; }

remora_dev=()
[[ "$DEVICE" != "cpu" ]] && remora_dev=(--device "$DEVICE")

# On any failure: name the step and show the end of its log (dorado/remora/modkit log to files).
STEP=""; STEP_LOG=""
on_err() {
  local rc=$?
  echo "ERROR (exit $rc) during: ${STEP:-unknown step}" >&2
  if [[ -n "$STEP_LOG" && -s "$STEP_LOG" ]]; then
    echo "--- last lines of $STEP_LOG ---" >&2
    tail -n 25 "$STEP_LOG" | tr '\r' '\n' | tail -n 25 >&2
  fi
  exit "$rc"
}
trap on_err ERR

mkdir -p "$OUT"
mapfile -t ROWS < <(tail -n +2 "$SAMPLES")

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r construct barcode _pos _offset _pod5 <<< "$row"
  pod5="$POD5_DIR/${construct}${POD5_SUFFIX}"
  [[ -f "$pod5" ]] || { echo "ERROR: missing $pod5" >&2; exit 1; }
  p="$OUT/$construct"
  echo "=== $construct ($barcode) ==="

  # 1. basecall with move table and inline alignment (same flags as the original run)
  STEP="$construct: dorado basecaller ($BASECALL_MODEL)"; STEP_LOG="$p.dorado.log"
  run_dorado basecaller "$BASECALL_MODEL" "$pod5" --models-directory "$MODELS_DIR" \
    --emit-moves --kit-name "$KIT" --reference "$REF" \
    > "$p.mapped.bam" 2> "$p.dorado.log"

  # 2. drop supplementary alignments (both strands kept), optionally restrict to the barcode
  STEP="$construct: samtools filter/sort/index"; STEP_LOG=""
  if [[ "$BC_FILTER" == "1" ]]; then
    samtools view -b -F 2048 -d "BC:${KIT}_${barcode}" "$p.mapped.bam" | samtools sort -o "$p.bam"
  else
    samtools view -b -F 2048 "$p.mapped.bam" | samtools sort -o "$p.bam"
  fi
  samtools index "$p.bam"
  n=$(samtools view -c "$p.bam")
  echo "  reads after filtering: $n"
  if [[ "$n" -eq 0 ]]; then
    echo "ERROR: 0 reads kept for $construct. If your dorado version writes a different BC tag," >&2
    echo "       rerun with BC_FILTER=0." >&2
    exit 1
  fi

  # 3. 8-oxoG inference, reference-anchored (this is how the published tables were produced)
  STEP="$construct: remora infer"; STEP_LOG="$p.remora.log"
  remora infer from_pod5_and_bam "$pod5" "$p.bam" \
    --model "$MODEL" --out-bam "$p.infer.bam" --reference-anchored ${remora_dev[@]+"${remora_dev[@]}"} \
    --log-filename "$p.remora.log"
  samtools sort -o "$p.infer.sorted.bam" "$p.infer.bam"
  samtools index "$p.infer.sorted.bam"

  # 4. pileups: unfiltered (window analysis, Table 4) and fixed-threshold (per-base, Table 3)
  STEP="$construct: modkit pileup --no-filtering"; STEP_LOG="$p.modkit_nofilt.log"
  modkit pileup "$p.infer.sorted.bam" "$p.nofilt.bed" --ref "$REF" --no-filtering \
    --log-filepath "$p.modkit_nofilt.log"
  STEP="$construct: modkit pileup --filter-threshold"; STEP_LOG="$p.modkit_filt.log"
  modkit pileup "$p.infer.sorted.bam" "$p.filt.bed" --ref "$REF" \
    --filter-threshold "G:${FILTER_THRESHOLD}" --log-filepath "$p.modkit_filt.log"
done

echo "Done. Pileups in $OUT (*.nofilt.bed, *.filt.bed)."
