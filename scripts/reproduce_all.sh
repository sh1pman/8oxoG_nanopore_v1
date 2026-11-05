#!/usr/bin/env bash
# One-shot reproduction: plasmid inference -> Tables 3 and 4 -> cassette ext-val.
#
# Usage:
#   bash scripts/reproduce_all.sh <pod5_dir> [out_dir]                 full data (Zenodo)
#   POD5_SUFFIX=_demo.pod5 bash scripts/reproduce_all.sh demo/pod5      quick demo (qualitative)

set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
POD5_DIR=${1:?usage: reproduce_all.sh <pod5_dir> [out_dir]}
OUT=${2:-results}

bash "$REPO/scripts/run_plasmid_inference.sh" "$POD5_DIR" "$OUT/plasmids"
python3 "$REPO/scripts/table3_per_base.py"   --results "$OUT/plasmids" --outdir "$OUT/tables"
python3 "$REPO/scripts/table4_window_scan.py" --results "$OUT/plasmids" --outdir "$OUT/tables"
bash "$REPO/scripts/validate_ext_val.sh" "$OUT/ext_val"

echo "All outputs in $OUT/. Compare $OUT/tables/*.tsv with expected/."
