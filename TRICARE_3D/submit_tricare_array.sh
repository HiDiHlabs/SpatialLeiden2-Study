#!/usr/bin/env bash
# Submit TRICARE tile (CPU) or univ2 (GPU) as a SLURM array: one specimen per task.
#
#   export TRICARE=./data/TRICARE
#   export COHORTS="prostate_test"
#   export MODEL_DIR=models          # univ2 only
#   bash TRICARE_3D/submit_tricare_array.sh tile
#   bash TRICARE_3D/submit_tricare_array.sh univ2 --skip-existing
#
# Optional: SLURM_PARTITION, SLURM_GPU_PARTITION, SLURM_TIME, SLURM_MEM, SLURM_CPUS,
#           CONDA_ENV, OUT_DIR, TRICARE_EXTRA_ARGS
set -euo pipefail

STAGE="${1:-}"
shift || true
# Remaining CLI flags go to the python worker (e.g. --skip-existing)
USER_EXTRA=("$@")

if [[ "$STAGE" != "tile" && "$STAGE" != "univ2" ]]; then
  echo "Usage: $0 tile|univ2 [extra python flags...]" >&2
  exit 2
fi

TRICARE="${TRICARE:?Set TRICARE (e.g. export TRICARE=./data/TRICARE)}"
COHORTS="${COHORTS:-prostate_test}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/run_tricare_univ2_multileiden.py"
OUT_DIR="${OUT_DIR:-${TRICARE}/univ2_multileiden}"
LOG_DIR="${LOG_DIR:-${OUT_DIR}/logs/slurm}"
mkdir -p "$LOG_DIR" "$OUT_DIR"

# Emit a bash snippet that defines conda for non-interactive jobs (no conda init needed).
_conda_bootstrap_snippet() {
  cat <<'SNIP'
_tricare_conda_activate() {
  local env="$1"
  local f
  for f in \
    "${HOME}/work/miniconda/etc/profile.d/conda.sh" \
    "${HOME}/miniconda3/etc/profile.d/conda.sh" \
    "${HOME}/miniforge3/etc/profile.d/conda.sh" \
    "${HOME}/mambaforge/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh"
  do
    if [[ -f "$f" ]]; then
      # shellcheck source=/dev/null
      source "$f"
      conda activate "$env"
      return 0
    fi
  done
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$env"
    return 0
  fi
  # Last resort: run env python directly if we can locate it next to conda
  local base bin
  base="$(command -v conda 2>/dev/null | xargs dirname 2>/dev/null | xargs dirname 2>/dev/null || true)"
  for base in "$base" "${HOME}/miniconda3" "${HOME}/miniforge3" "${HOME}/mambaforge" "${HOME}/work/miniconda"; do
    bin="${base}/envs/${env}/bin"
    if [[ -x "${bin}/python" ]]; then
      export PATH="${bin}:${PATH}"
      return 0
    fi
  done
  echo "Could not activate conda env: $env" >&2
  return 1
}
SNIP
}

_resolve_univ2_model_dir() {
  # Prefer explicit MODEL_DIR / UNIV2_MODEL_DIR; else search common locations.
  local cand abs
  local -a cands=()
  if [[ -n "${MODEL_DIR:-}" ]]; then
    cands+=("$MODEL_DIR")
  fi
  if [[ -n "${UNIV2_MODEL_DIR:-}" ]]; then
    cands+=("$UNIV2_MODEL_DIR")
  fi
  cands+=(
    "${REPO_ROOT}/models"
    "${REPO_ROOT}/models/uni2-h"
    "${PWD}/models"
    "${PWD}/models/uni2-h"
    "${HOME}/models"
    "${HOME}/models/uni2-h"
    "${HOME}/multileiden/models"
    "${HOME}/multileiden/models/uni2-h"
    "${TRICARE}/../models"
    "${TRICARE}/../models/uni2-h"
    "${TRICARE}/models"
  )
  for cand in "${cands[@]}"; do
    [[ -z "$cand" ]] && continue
    if [[ "$cand" = /* ]]; then
      abs="$cand"
    else
      abs="$(cd "$REPO_ROOT" 2>/dev/null && realpath -m "$cand" 2>/dev/null || echo "${REPO_ROOT}/${cand}")"
    fi
    if [[ -f "${abs}/pytorch_model.bin" ]]; then
      echo "$abs"
      return 0
    fi
  done
  # Hugging Face hub cache (timm / hf_hub_download)
  local snap
  snap="$(find "${HOME}/.cache/huggingface/hub/models--MahmoodLab--UNI2-h" -name pytorch_model.bin 2>/dev/null | head -1 || true)"
  if [[ -n "$snap" && -f "$snap" ]]; then
    echo "$(dirname "$snap")"
    return 0
  fi
  return 1
}

if [[ "$STAGE" == "tile" ]]; then
  CONDA_ENV_NAME="${CONDA_ENV:-tricare_tile}"
  PARTITION="${SLURM_PARTITION:-standard}"
  TIME="${SLURM_TIME:-12:00:00}"
  MEM="${SLURM_MEM:-32G}"
  CPUS="${SLURM_CPUS:-4}"
  GRES_LINE=""
  STAGE_ARGS=(--tile-backend lazyslide --find-tissues --skip-existing)
else
  CONDA_ENV_NAME="${CONDA_ENV:-tricare_univ2}"
  PARTITION="${SLURM_GPU_PARTITION:-gpu}"
  TIME="${SLURM_TIME:-12:00:00}"
  MEM="${SLURM_MEM:-32G}"
  CPUS="${SLURM_CPUS:-4}"
  GRES_LINE="#SBATCH --gres=gpu:1"
  if ! MODEL_DIR="$(_resolve_univ2_model_dir)"; then
    echo "UNI2 weights not found (need DIRECTORY/pytorch_model.bin)." >&2
    echo "They are gated on Hugging Face — download once, then resubmit:" >&2
    echo "  # request access: https://huggingface.co/MahmoodLab/UNI2-h" >&2
    echo "  conda activate tricare_univ2" >&2
    echo "  python -m huggingface_hub.commands.huggingface_cli login" >&2
    echo "  export MODEL_DIR=./models" >&2
    echo "  python TRICARE_3D/download_univ2_weights.py --model-dir \"\$MODEL_DIR\"" >&2
    echo "  bash TRICARE_3D/submit_tricare_array.sh univ2" >&2
    exit 1
  fi
  echo "UNI2 model-dir: $MODEL_DIR"
  STAGE_ARGS=(--model-dir "$MODEL_DIR" --device auto --skip-existing)
fi

# Submit-time listing: use current python (caller often already has conda env active).
# Do not call `conda activate` here — non-interactive bash lacks the conda shell hook.
if ! python -c "import mudata" 2>/dev/null; then
  echo "Need python with mudata on PATH. Activate tricare_tile or download_tricare, then re-run." >&2
  exit 1
fi
# shellcheck disable=SC2206
COHORT_ARR=($COHORTS)
LIST="${OUT_DIR}/specimen_list.txt"
python "$PY_SCRIPT" \
  --tricare-dir "$TRICARE" \
  --out-dir "$OUT_DIR" \
  --cohorts "${COHORT_ARR[@]}" \
  --list-specimens > "$LIST"
N="$(wc -l < "$LIST" | tr -d ' ')"
if [[ "$N" -lt 1 ]]; then
  echo "No specimens listed in $LIST" >&2
  exit 1
fi
MAX_IDX=$((N - 1))
echo "Submitting $STAGE array 0-${MAX_IDX} ($N specimens) partition=$PARTITION → $LOG_DIR"

ARGS_FILE="${LOG_DIR}/${STAGE}_args.txt"
{
  printf '%s\0' --tricare-dir "$TRICARE"
  printf '%s\0' --out-dir "$OUT_DIR"
  printf '%s\0' --cohorts
  # shellcheck disable=SC2086
  for c in $COHORTS; do printf '%s\0' "$c"; done
  printf '%s\0' --stage "$STAGE"
  printf '%s\0' "${STAGE_ARGS[@]}"
  if ((${#USER_EXTRA[@]})); then
    printf '%s\0' "${USER_EXTRA[@]}"
  fi
  if [[ -n "${TRICARE_EXTRA_ARGS:-}" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARR=($TRICARE_EXTRA_ARGS)
    printf '%s\0' "${EXTRA_ARR[@]}"
  fi
} > "$ARGS_FILE"

JOB_SCRIPT="${LOG_DIR}/_worker_${STAGE}.sh"
{
  cat <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=tricare_${STAGE}
#SBATCH --partition=${PARTITION}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${TIME}
#SBATCH --array=0-${MAX_IDX}
#SBATCH --output=${LOG_DIR}/${STAGE}_%A_%a.out
#SBATCH --error=${LOG_DIR}/${STAGE}_%A_%a.err
${GRES_LINE}
set -euo pipefail
EOF
  _conda_bootstrap_snippet
  cat <<EOF
_tricare_conda_activate ${CONDA_ENV_NAME}
cd "${REPO_ROOT}"
IDX="\${SLURM_ARRAY_TASK_ID}"
mapfile -d '' -t ARGS < "${ARGS_FILE}"
echo "[tricare ${STAGE}] specimen-index=\${IDX} host=\$(hostname) python=\$(command -v python)"
python "${PY_SCRIPT}" "\${ARGS[@]}" --specimen-index "\${IDX}"
EOF
} > "$JOB_SCRIPT"

chmod +x "$JOB_SCRIPT"
sbatch "$JOB_SCRIPT"
echo "Specimen list: $LIST"
echo "Worker script: $JOB_SCRIPT"
