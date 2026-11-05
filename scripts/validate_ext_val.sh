#!/usr/bin/env bash
# Reproduce the per-chunk validation of model v1 on the held-out cassette sets
# (canonical G vs 8-oxoG at the same site and anchor context; manuscript: 99.9%).
#
# Usage:  bash scripts/validate_ext_val.sh [out_dir]          (default: results/ext_val)
#   DEVICE=0 (GPU id) or DEVICE=cpu
#
# Remora writes absolute paths into dataset configs, so the config is built here, on your machine.

set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-results/ext_val}
MODEL="$REPO/model/remora_8oxoG_v1.pt"
NEG="$REPO/ext_val/ext_val_neg"     # Cas_G library, site 44, canonical G
POS="$REPO/ext_val/ext_val_pos"     # Cas_oG library, site 44, 8-oxoG
DEVICE=${DEVICE:-0}
dev=(); [[ "$DEVICE" != "cpu" ]] && dev=(--device "$DEVICE")

mkdir -p "$OUT"

# 1. Balanced set, exactly as evaluated during training: batches of 2048 chunks drawn 1:1 from both
#    sets until the smaller set cannot fill a batch -> 1024 + 1024 chunks.
remora dataset make_config "$OUT/ext_val.jsn" "$NEG" "$POS" --dataset-weights 1 1 \
  --log-filename "$OUT/make_config.log"
remora validate from_remora_dataset "$OUT/ext_val.jsn" --model "$MODEL" \
  --batch-size 2048 --out-file "$OUT/ext_val_balanced.txt" \
  --log-filename "$OUT/validate_balanced.log" ${dev[@]+"${dev[@]}"}

#    Note: remora splits each batch between the two sets at random (multinomial), so the class counts
#    vary between runs (e.g. 1022/1026 instead of 1024/1024); use step 2 for deterministic numbers.

# 2. Every chunk of each set separately (specificity on canonical G, sensitivity on 8-oxoG).
#    Remora drops an incomplete final batch, so batch size 1 is used to score every chunk.
for ds in ext_val_neg ext_val_pos; do
  remora validate from_remora_dataset "$REPO/ext_val/$ds" --model "$MODEL" \
    --batch-size 1 --out-file "$OUT/${ds}_all.txt" --log-filename "$OUT/validate_${ds}.log" ${dev[@]+"${dev[@]}"}
done

echo "Results:"
for f in "$OUT"/ext_val_balanced.txt "$OUT"/ext_val_neg_all.txt "$OUT"/ext_val_pos_all.txt; do
  echo "--- $f"; cat "$f"
done
