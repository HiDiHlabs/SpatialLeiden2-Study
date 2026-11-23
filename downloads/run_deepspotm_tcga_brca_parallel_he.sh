#!/usr/bin/env bash
# Disk-friendly DeepSpot-M TCGA-BRCA: select → PCC → per-sample
#   (atlas h5ad ↓ → SVS ↓ → HE extract) → MuData. SVS kept by default
#   (--delete-svs to remove after HE).
#
# Prerequisite: envs download_deepspotm_tcga + lazyslide + univ2; huggingface-cli login;
# Prefer: bash downloads/submit_deepspotm_tcga_brca_slurm.sh he
#               UNI2: MODEL_DIR directory with pytorch_model.bin (default: ./models).
#
# From repo root:
#   bash downloads/run_deepspotm_tcga_brca_parallel_he.sh --help
#   bash downloads/run_deepspotm_tcga_brca_parallel_he.sh \
#     --out-dir ./data/TCGA_BRCA --max-patients 100 --jobs 1
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# ── defaults (override via flags or env) ─────────────────────────────────────
OUT_DIR="./data/TCGA_BRCA"
MAX_PATIENTS=100
MAX_SLIDES=""
MIN_SPOTS=500
STAGES="II,III"
ENRICH_PFI="1:1"
ENRICH_SEED=0
REQUIRE_ENDPOINT="PFI"
JOBS=1                        # disk-safe: one SVS at a time
HE_ENCODER="univ2"
UNIV2_MODEL_DIR_ARG=""        # set via --univ2-model-dir; else ./models
HE_DEVICE="${HE_DEVICE:-auto}"
HE_BATCH_SIZE=32
HE_NPY_DIR=""
SVS_DIR=""
DELETE_SVS=0
BULK_ATLAS=0                  # 0=fetch h5ad per sample in HE; 1=prefetch all
RUN_PCC=0                     # 0=keep all genes (default); 1=run PCC filter
PCC_CAMEO=1                   # only used when RUN_PCC=1
PCC_THRESHOLD=0.4
PCC_DEVICE="cpu"
PCC_MAX_SPOTS=2000
PCC_CAMEO_GENES_CSV=""
GENE_PCC_CSV=""
PCC_VISIUM_H5AD=""
PCC_TILE_DIR=""
SKIP_SELECT=0
SKIP_HE=0
SKIP_MUDATA=0
RESELECT=0
PY="${PYTHON:-python}"

# Default UNI2 weights dir: ./models (override with MODEL_DIR / --univ2-model-dir)
_resolve_univ2_model_dir() {
  local cand abs
  local -a cands=()
  [[ -n "${UNIV2_MODEL_DIR_ARG}" ]] && cands+=("$UNIV2_MODEL_DIR_ARG")
  [[ -n "${MODEL_DIR:-}" ]] && cands+=("$MODEL_DIR")
  [[ -n "${UNIV2_MODEL_DIR:-}" ]] && cands+=("$UNIV2_MODEL_DIR")
  cands+=("${ROOT}/models")
  for cand in "${cands[@]}"; do
    [[ -z "$cand" ]] && continue
    if [[ "$cand" = /* ]]; then
      abs="$cand"
    else
      abs="$(cd "$ROOT" 2>/dev/null && realpath -m "$cand" 2>/dev/null || echo "${ROOT}/${cand}")"
    fi
    if [[ -f "${abs}/pytorch_model.bin" ]]; then
      echo "$abs"
      return 0
    fi
  done
  return 1
}

usage() {
  cat <<'EOF'
Disk-friendly SVS → HE → delete SVS for DeepSpot-M TCGA-BRCA (then MuData).

Per sample (default --jobs 1): download atlas h5ad (if needed) → download SVS
→ UNI2 HE → delete SVS → next sample. No bulk SVS queue.

Usage:
  bash downloads/run_deepspotm_tcga_brca_parallel_he.sh [flags]

I/O
  --out-dir DIR              Output root (default: ./data/TCGA_BRCA)
  --he-npy-dir DIR           HE .npy dir (default: {out}/he_embeddings)
  --svs-dir DIR              SVS dir (default: {out}/svs)

Clinical select
  --stages LIST              Collapsed AJCC stages (default: II,III; or all)
  --enrich-pfi E:C           PFI event:censored ratio (default: 1:1)
  --no-enrich-pfi            Disable PFI enrichment
  --enrich-seed N            (default: 0)
  --max-patients N           (default: 100)
  --max-slides N
  --min-spots N              (default: 500)
  --require-endpoint NAME    any|OS|PFI|both (default: PFI)
  --reselect

HE
  --jobs|-j N                Concurrent samples (default: 1)
  --he-encoder NAME          univ2|deepspotm_midnight|mean_rgb (default: univ2)
  --univ2-model-dir DIR      UNI2 dir with pytorch_model.bin
                             (default: ./models)
  --he-device NAME            auto|cuda|cpu
  --he-batch-size N
  --keep-svs                 Keep SVS after HE (default)
  --delete-svs               Delete SVS after HE
  --bulk-atlas               Prefetch all selected atlas h5ads before HE
  --no-bulk-atlas            Per-sample atlas fetch during HE (default)

PCC (off by default — MuData keeps all genes)
  --pcc-filter               Enable PCC gene subsetting
  --skip-pcc                 Disable PCC (default)
  --pcc-cameo / --no-pcc-cameo
  --pcc-threshold X
  --pcc-device NAME
  --pcc-max-spots N
  --pcc-cameo-genes-csv PATH
  --gene-pcc-csv PATH
  --pcc-visium-h5ad PATH
  --pcc-tile-dir PATH

Stage skips
  --skip-select|--skip-he|--skip-mudata

Other
  --python PATH
  -h|--help

Pipeline: select → [optional PCC] → per-sample (h5ad↓ → SVS↓ → HE) → mudata
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --he-npy-dir) HE_NPY_DIR="$2"; shift 2 ;;
    --svs-dir) SVS_DIR="$2"; shift 2 ;;
    --stages) STAGES="$2"; shift 2 ;;
    --enrich-pfi) ENRICH_PFI="$2"; shift 2 ;;
    --no-enrich-pfi) ENRICH_PFI=""; shift ;;
    --enrich-seed) ENRICH_SEED="$2"; shift 2 ;;
    --max-patients) MAX_PATIENTS="$2"; shift 2 ;;
    --max-slides) MAX_SLIDES="$2"; shift 2 ;;
    --min-spots) MIN_SPOTS="$2"; shift 2 ;;
    --require-endpoint) REQUIRE_ENDPOINT="$2"; shift 2 ;;
    --reselect) RESELECT=1; shift ;;
    --jobs|-j) JOBS="$2"; shift 2 ;;
    --he-encoder) HE_ENCODER="$2"; shift 2 ;;
    --univ2-model-dir) UNIV2_MODEL_DIR_ARG="$2"; shift 2 ;;
    --he-device) HE_DEVICE="$2"; shift 2 ;;
    --he-batch-size) HE_BATCH_SIZE="$2"; shift 2 ;;
    --delete-svs) DELETE_SVS=1; shift ;;
    --keep-svs) DELETE_SVS=0; shift ;;
    --bulk-atlas) BULK_ATLAS=1; shift ;;
    --no-bulk-atlas) BULK_ATLAS=0; shift ;;
    --pcc-filter) RUN_PCC=1; shift ;;
    --skip-pcc) RUN_PCC=0; shift ;;
    --pcc-cameo) PCC_CAMEO=1; RUN_PCC=1; shift ;;
    --no-pcc-cameo) PCC_CAMEO=0; shift ;;
    --pcc-threshold) PCC_THRESHOLD="$2"; shift 2 ;;
    --pcc-device) PCC_DEVICE="$2"; shift 2 ;;
    --pcc-max-spots) PCC_MAX_SPOTS="$2"; shift 2 ;;
    --pcc-cameo-genes-csv) PCC_CAMEO_GENES_CSV="$2"; shift 2 ;;
    --gene-pcc-csv) GENE_PCC_CSV="$2"; RUN_PCC=1; shift 2 ;;
    --pcc-visium-h5ad) PCC_VISIUM_H5AD="$2"; RUN_PCC=1; shift 2 ;;
    --pcc-tile-dir) PCC_TILE_DIR="$2"; shift 2 ;;
    --skip-select) SKIP_SELECT=1; shift ;;
    --skip-atlas) BULK_ATLAS=0; shift ;;  # alias: never bulk-prefetch
    --skip-he) SKIP_HE=1; shift ;;
    --skip-mudata) SKIP_MUDATA=1; shift ;;
    --python) PY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1 (try --help)" >&2; exit 1 ;;
  esac
done

SCRIPT="downloads/retrieve_deepspotm_tcga_brca.py"
SEL="${OUT_DIR}/selected_slides.csv"
IDS="${OUT_DIR}/selected_sample_ids.txt"
LOG_DIR="${OUT_DIR}/logs/he"
mkdir -p "$OUT_DIR" "$LOG_DIR"

if ! UNIV2_MODEL_DIR="$(_resolve_univ2_model_dir)"; then
  echo "UNI2 weights not found (need DIRECTORY/pytorch_model.bin)." >&2
  echo "Default: ${ROOT}/models (or set MODEL_DIR / --univ2-model-dir)." >&2
  echo "  ls ${ROOT}/models/pytorch_model.bin" >&2
  echo "Or: --univ2-model-dir /path/to/dir_with_pytorch_model.bin" >&2
  exit 1
fi
echo "UNI2 model-dir: $UNIV2_MODEL_DIR"

COMMON=(
  --out-dir "$OUT_DIR"
  --stages "$STAGES"
  --min-spots "$MIN_SPOTS"
  --he-encoder "$HE_ENCODER"
  --univ2-model-dir "$UNIV2_MODEL_DIR"
  --he-device "$HE_DEVICE"
  --he-batch-size "$HE_BATCH_SIZE"
)
[[ -n "$HE_NPY_DIR" ]] && COMMON+=(--he-npy-dir "$HE_NPY_DIR")
[[ -n "$SVS_DIR" ]] && COMMON+=(--svs-dir "$SVS_DIR")

SELECT_ARGS=(
  --stage select
  --require-endpoint "$REQUIRE_ENDPOINT"
  --enrich-seed "$ENRICH_SEED"
)
[[ -n "$MAX_PATIENTS" ]] && SELECT_ARGS+=(--max-patients "$MAX_PATIENTS")
[[ -n "$MAX_SLIDES" ]] && SELECT_ARGS+=(--max-slides "$MAX_SLIDES")
[[ -n "$ENRICH_PFI" ]] && SELECT_ARGS+=(--enrich-pfi "$ENRICH_PFI")
[[ "$RESELECT" -eq 1 ]] && SELECT_ARGS+=(--reselect)

PCC_ARGS=(--stage pcc --pcc-filter --pcc-threshold "$PCC_THRESHOLD" --pcc-device "$PCC_DEVICE" --pcc-max-spots "$PCC_MAX_SPOTS")
[[ "$PCC_CAMEO" -eq 1 ]] && PCC_ARGS+=(--pcc-cameo)
[[ -n "$PCC_CAMEO_GENES_CSV" ]] && PCC_ARGS+=(--pcc-cameo-genes-csv "$PCC_CAMEO_GENES_CSV")
[[ -n "$GENE_PCC_CSV" ]] && PCC_ARGS+=(--gene-pcc-csv "$GENE_PCC_CSV")
[[ -n "$PCC_VISIUM_H5AD" ]] && PCC_ARGS+=(--pcc-visium-h5ad "$PCC_VISIUM_H5AD")
[[ -n "$PCC_TILE_DIR" ]] && PCC_ARGS+=(--pcc-tile-dir "$PCC_TILE_DIR")

echo "Config: out=${OUT_DIR} max-patients=${MAX_PATIENTS} jobs=${JOBS} he=${HE_ENCODER} delete-svs=${DELETE_SVS} pcc=${RUN_PCC}"

if [[ "$SKIP_SELECT" -eq 0 ]]; then
  echo "[1/4] select"
  "$PY" "$SCRIPT" "${COMMON[@]}" "${SELECT_ARGS[@]}"
else
  echo "[1/4] skip select (reuse ${SEL})"
fi

if [[ ! -f "$SEL" ]]; then
  echo "Missing ${SEL}; run without --skip-select" >&2
  exit 1
fi

"$PY" - <<PY
import pandas as pd
from pathlib import Path
df = pd.read_csv("${SEL}")
Path("${IDS}").write_text("\n".join(df["sample_id"].astype(str)) + "\n")
print(f"Wrote ${IDS} ({len(df)} samples)")
PY

if [[ "$RUN_PCC" -eq 1 ]]; then
  echo "[2/4] PCC gene filter"
  "$PY" "$SCRIPT" "${COMMON[@]}" "${PCC_ARGS[@]}"
else
  echo "[2/4] skip PCC (keep all genes)"
fi

if [[ "$BULK_ATLAS" -eq 1 ]]; then
  echo "[2b] bulk atlas prefetch (all selected h5ads)"
  "$PY" "$SCRIPT" "${COMMON[@]}" --stage download
fi

he_one() {
  local sid="$1"
  local log="${LOG_DIR}/${sid}.log"
  local npy_dir="${HE_NPY_DIR:-${OUT_DIR}/he_embeddings}"
  local svs_dir="${SVS_DIR:-${OUT_DIR}/svs}"
  echo "[HE] ${sid}"
  if [[ -f "${npy_dir}/${sid}.npy" ]]; then
    echo "[HE] skip (npy exists): ${sid}"
    if [[ "${DELETE_SVS}" -eq 1 ]]; then
      rm -f "${svs_dir}/${sid}.svs"
    fi
    return 0
  fi
  local cmd=(
    "$PY" "$SCRIPT"
    --out-dir "$OUT_DIR"
    --stage he
    --sample-id "$sid"
    --download-atlas
    --download-svs
    --he-encoder "$HE_ENCODER"
    --univ2-model-dir "$UNIV2_MODEL_DIR"
    --he-device "$HE_DEVICE"
    --he-batch-size "$HE_BATCH_SIZE"
  )
  [[ -n "$HE_NPY_DIR" ]] && cmd+=(--he-npy-dir "$HE_NPY_DIR")
  [[ -n "$SVS_DIR" ]] && cmd+=(--svs-dir "$SVS_DIR")
  [[ "$DELETE_SVS" -eq 1 ]] && cmd+=(--delete-svs)
  set +e
  "${cmd[@]}" >"$log" 2>&1
  local rc=$?
  set -e
  if [[ "$rc" -ne 0 ]]; then
    echo "[HE] FAILED ${sid} (exit ${rc}). Last log lines:" >&2
    tail -n 40 "$log" >&2 || true
    echo "Full log: ${log}" >&2
    return "$rc"
  fi
  if [[ "${DELETE_SVS}" -eq 1 ]]; then
    rm -f "${svs_dir}/${sid}.svs"
  fi
}

export -f he_one
export OUT_DIR SCRIPT PY HE_ENCODER UNIV2_MODEL_DIR HE_DEVICE HE_BATCH_SIZE
export HE_NPY_DIR SVS_DIR DELETE_SVS LOG_DIR ROOT PATH

if [[ "$SKIP_HE" -eq 0 ]]; then
  echo "[3/4] HE per sample: h5ad↓ → SVS↓ → ${HE_ENCODER} (jobs=${JOBS}; delete-svs=${DELETE_SVS})"
  if [[ "$JOBS" -le 1 ]]; then
    while IFS= read -r sid || [[ -n "$sid" ]]; do
      [[ -z "$sid" ]] && continue
      he_one "$sid"
    done < "$IDS"
  elif command -v parallel >/dev/null 2>&1; then
    parallel -j "$JOBS" --halt soon,fail=1 he_one :::: "$IDS"
  else
    echo "GNU parallel not found; using xargs -P ${JOBS}"
    xargs -P "$JOBS" -n 1 -a "$IDS" bash -c 'he_one "$1"' _
  fi
else
  echo "[3/4] skip HE"
fi

if [[ "$SKIP_MUDATA" -eq 0 ]]; then
  echo "[4/4] assemble MuData"
  "$PY" "$SCRIPT" "${COMMON[@]}" --stage mudata
else
  echo "[4/4] skip mudata"
fi

echo "Done. HE: ${HE_NPY_DIR:-${OUT_DIR}/he_embeddings}/  MuData: ${OUT_DIR}/mudata/"
svs_root="${SVS_DIR:-${OUT_DIR}/svs}"
leftover=$(find "$svs_root" -name '*.svs' 2>/dev/null | wc -l | tr -d ' ')
echo "Leftover SVS files: ${leftover}"
