#!/usr/bin/env bash
# SLURM submit for DeepSpot-M TCGA-BRCA: select (CPU) → HE array (GPU) → MuData (CPU).
# Each HE array task uses project conda envs:
#   fetch  → download_deepspotm_tcga (HF atlas + GDC SVS)
#   tile   → lazyslide (non_overlap crop by default; CROP_SIZE=hest|…)
#   univ2  → univ2 (univ2_snake.py)
# then writes he_embeddings/{id}.npy + he_overview/{id}.jpg; keeps SVS; deletes tiles by default.
#
# Array task logs: ${OUT_DIR}/logs/slurm/<array_job_id>/he_<task>.{out,err}
# Mudata (when chained): same directory as the HE array job.
#
# From repo root (after HF login + UNI2 weights):
#   export OUT_DIR=./data/TCGA_BRCA
#   export MODEL_DIR=models
#   bash downloads/submit_deepspotm_tcga_brca_slurm.sh all
#
# Stages:
#   all | select | he | mudata
#   test [SAMPLE_ID]   # fetch→tile→univ2 locally for one sample (debug; keeps tiles)
#
# Optional env:
#   OUT_DIR, MAX_PATIENTS, ENRICH_PFI, STAGES
#   MODEL_DIR / UNIV2_MODEL_DIR (default: ./models)
#   CONDA_ENV (select/mudata/fetch; default: download_deepspotm_tcga)
#   CONDA_ENV_TILE (default: lazyslide), CONDA_ENV_HE (default: univ2)
#   SLURM_GPU_PARTITION (default: gpu), SLURM_PARTITION (CPU; default: standard)
#   SLURM_TIME, SLURM_MEM, SLURM_CPUS, SLURM_HE_TIME, SLURM_HE_MEM
#   SLURM_ARRAY_THROTTLE (optional; e.g. 4 → --array=0-N%4; unset = no cap)
#   DELETE_TILES=0 to keep JPEG tiles after UNI2 (test defaults to keep)
#   DELETE_SVS=1 to remove SVS after HE (default: keep)
#   FORCE_HE=1 to redo UNI2 even if npy exists (test defaults to 1)
#   FORCE_RETILE=1 to re-run lazyslide even if tiles/{id}.h5ad exists
#   CROP_SIZE (default: non_overlap — Visium √3/2 pitch; hest=112 µm overlaps)
#   RESELECT=1, HE_BATCH_SIZE
set -euo pipefail

STAGE="${1:-all}"
shift || true
USER_EXTRA=("$@")

case "$STAGE" in
  all|select|he|mudata|test) ;;
  -h|--help)
    sed -n '2,36p' "$0"
    exit 0
    ;;
  *)
    echo "Usage: $0 all|select|he|mudata|test [SAMPLE_ID]" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$REPO_ROOT"

OUT_DIR="${OUT_DIR:-./data/TCGA_BRCA}"
MAX_PATIENTS="${MAX_PATIENTS:-100}"
ENRICH_PFI="${ENRICH_PFI:-1:1}"
STAGES_AJCC="${STAGES:-II,III}"
CONDA_ENV_NAME="${CONDA_ENV:-download_deepspotm_tcga}"
CONDA_ENV_TILE_NAME="${CONDA_ENV_TILE:-lazyslide}"
CONDA_ENV_HE_NAME="${CONDA_ENV_HE:-univ2}"
CROP_SIZE="${CROP_SIZE:-non_overlap}"
LOG_ROOT="${LOG_DIR:-${OUT_DIR}/logs/slurm}"
IDS="${OUT_DIR}/selected_sample_ids.txt"
SEL="${OUT_DIR}/selected_slides.csv"
PY_SCRIPT="${SCRIPT_DIR}/retrieve_deepspotm_tcga_brca.py"
TILE_SCRIPT="${SCRIPT_DIR}/lazyslide_tile.py"
UNIV2_SCRIPT="${SCRIPT_DIR}/univ2_snake.py"
mkdir -p "$LOG_ROOT" "$OUT_DIR"

CPU_PARTITION="${SLURM_PARTITION:-standard}"
GPU_PARTITION="${SLURM_GPU_PARTITION:-gpu}"
CPU_TIME="${SLURM_TIME:-4:00:00}"
HE_TIME="${SLURM_HE_TIME:-12:00:00}"
CPU_MEM="${SLURM_MEM:-16G}"
HE_MEM="${SLURM_HE_MEM:-32G}"
CPUS="${SLURM_CPUS:-4}"
HE_CPUS="${SLURM_HE_CPUS:-4}"
ARRAY_THROTTLE="${SLURM_ARRAY_THROTTLE:-}"

_deepspotm_conda_activate() {
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

_conda_bootstrap_snippet() {
  # Emit the activate function into SLURM worker scripts.
  declare -f _deepspotm_conda_activate
}

_resolve_univ2_model_dir() {
  local cand abs
  local -a cands=()
  [[ -n "${MODEL_DIR:-}" ]] && cands+=("$MODEL_DIR")
  [[ -n "${UNIV2_MODEL_DIR:-}" ]] && cands+=("$UNIV2_MODEL_DIR")
  cands+=("${REPO_ROOT}/models")
  for cand in "${cands[@]}"; do
    [[ -z "$cand" ]] && continue
    if [[ "$cand" = /* ]]; then abs="$cand"; else abs="${REPO_ROOT}/${cand}"; fi
    if [[ -f "${abs}/pytorch_model.bin" ]]; then
      echo "$abs"
      return 0
    fi
  done
  return 1
}

if ! MODEL_DIR_RESOLVED="$(_resolve_univ2_model_dir)"; then
  echo "UNI2 weights not found (need DIRECTORY/pytorch_model.bin)." >&2
  echo "Default search: ./models next to the study repo root, or set MODEL_DIR." >&2
  echo "  ls ${REPO_ROOT}/models/pytorch_model.bin" >&2
  echo "Or: export MODEL_DIR=/path/to/dir_with_pytorch_model.bin" >&2
  exit 1
fi
echo "UNI2 model-dir: $MODEL_DIR_RESOLVED"

_run_select_local() {
  local py="${PYTHON:-python}"
  if ! "$py" -c "import pandas" 2>/dev/null; then
    _deepspotm_conda_activate "$CONDA_ENV_NAME"
    py=python
  fi
  local sel_args=(
    --out-dir "$OUT_DIR"
    --stage select
    --stages "$STAGES_AJCC"
    --require-endpoint PFI
    --max-patients "$MAX_PATIENTS"
    --univ2-model-dir "$MODEL_DIR_RESOLVED"
  )
  [[ -n "$ENRICH_PFI" ]] && sel_args+=(--enrich-pfi "$ENRICH_PFI")
  [[ "${RESELECT:-0}" == "1" ]] && sel_args+=(--reselect)
  if ((${#USER_EXTRA[@]})); then
    sel_args+=("${USER_EXTRA[@]}")
  fi
  echo "[select] writing ${SEL}"
  "$py" "$PY_SCRIPT" "${sel_args[@]}"
  "$py" - <<PY
import pandas as pd
from pathlib import Path
df = pd.read_csv("${SEL}")
Path("${IDS}").write_text("\\n".join(df["sample_id"].astype(str)) + "\\n")
print(f"Wrote ${IDS} ({len(df)} samples)")
PY
}

# Run fetch → tile → univ2 for one sample (SID).
# Reuses tiles/{sid}.h5ad when present unless FORCE_RETILE=1.
_save_he_overview_if_needed() {
  local sid="$1"
  local svs="${OUT_DIR}/svs/${sid}.svs"
  local overview="${OUT_DIR}/he_overview/${sid}.jpg"
  if [[ -f "$overview" ]]; then
    return 0
  fi
  if [[ ! -f "$svs" ]]; then
    echo "[HE] no SVS for overview: ${sid}" >&2
    return 0
  fi
  _deepspotm_conda_activate "$CONDA_ENV_NAME"
  python - <<PY
from pathlib import Path
import sys
sys.path.insert(0, "${SCRIPT_DIR}")
from retrieve_deepspotm_tcga_brca import save_he_overview
save_he_overview(Path("${svs}"), Path("${overview}"))
PY
}

_run_he_one_sample() {
  local sid="$1"
  local npy="${OUT_DIR}/he_embeddings/${sid}.npy"
  local svs="${OUT_DIR}/svs/${sid}.svs"
  local tiles_dir="${OUT_DIR}/tiles/${sid}"
  local tiles_h5ad="${OUT_DIR}/tiles/${sid}.h5ad"
  local mf_h5ad="${OUT_DIR}/features/univ2/${sid}_mf.h5ad"
  local rel atlas_h5ad
  local need_tile=1

  if [[ -f "$npy" && "${FORCE_HE:-0}" != "1" ]]; then
    echo "[HE] skip (npy exists): ${sid}  (FORCE_HE=1 to redo UNI2)"
    _save_he_overview_if_needed "$sid"
    if [[ "${DELETE_SVS:-0}" == "1" ]]; then
      rm -f "$svs"
    fi
    return 0
  fi
  echo "[HE] sample=${sid} host=$(hostname) cuda=${CUDA_VISIBLE_DEVICES:-none}"

  if [[ -f "$tiles_h5ad" && "${FORCE_RETILE:-0}" != "1" ]]; then
    need_tile=0
    echo "[HE] reuse existing tiles: ${tiles_h5ad}  (FORCE_RETILE=1 to retile)"
  fi

  if [[ "$need_tile" -eq 1 ]]; then
    _deepspotm_conda_activate "$CONDA_ENV_NAME"
    python "$PY_SCRIPT" \
      --out-dir "$OUT_DIR" \
      --stage fetch \
      --sample-id "$sid"
    rel="$(python - <<PY
import pandas as pd
df = pd.read_csv("${SEL}")
row = df.loc[df["sample_id"].astype(str) == "${sid}"].iloc[0]
print(row["file_path"])
PY
)"
    atlas_h5ad="${OUT_DIR}/atlas/${rel}"
    if [[ ! -f "$atlas_h5ad" ]]; then
      echo "Missing atlas h5ad: ${atlas_h5ad}" >&2
      return 1
    fi
    if [[ ! -f "$svs" ]]; then
      echo "Missing SVS after fetch: ${svs}" >&2
      return 1
    fi

    _save_he_overview_if_needed "$sid"

    _deepspotm_conda_activate "$CONDA_ENV_TILE_NAME"
    python "$TILE_SCRIPT" \
      -s "$sid" \
      --wsi "$svs" \
      --coords-h5ad "$atlas_h5ad" \
      --crop_size "${CROP_SIZE}" \
      --tiles "$tiles_dir" \
      --adata_out "$tiles_h5ad"
  else
    _save_he_overview_if_needed "$sid"
  fi

  if [[ ! -f "$tiles_h5ad" ]]; then
    echo "Missing tiles h5ad: ${tiles_h5ad}" >&2
    return 1
  fi

  _deepspotm_conda_activate "$CONDA_ENV_HE_NAME"
  python - <<'PY'
import os, socket, torch
host = socket.gethostname()
built = getattr(torch.version, "cuda", None)
ok = torch.cuda.is_available()
print(f"torch={torch.__version__} cuda_built={built} cuda_available={ok} host={host} CVD={os.environ.get('CUDA_VISIBLE_DEVICES', '')!r}")
if not ok:
    raise SystemExit(
        "No CUDA device for UNI2. Run on a GPU node (sbatch/srun --gres=gpu:1), "
        f"not a CPU host like {host}. If on GPU and still failing, univ2 env may be CPU-only torch."
    )
PY
  mkdir -p "${OUT_DIR}/features/univ2"
  python "$UNIV2_SCRIPT" \
    --input "$tiles_h5ad" \
    --output "$mf_h5ad" \
    --device cuda \
    --batch-size "${HE_BATCH_SIZE:-32}" \
    --model-dir "$MODEL_DIR_RESOLVED"
  mkdir -p "${OUT_DIR}/he_embeddings"
  python - <<PY
import numpy as np
import scanpy as sc
from pathlib import Path
adata = sc.read_h5ad("${mf_h5ad}")
key = "univ2" if "univ2" in adata.obsm else next(iter(adata.obsm.keys()))
emb = np.asarray(adata.obsm[key], dtype=np.float32)
Path("${npy}").parent.mkdir(parents=True, exist_ok=True)
np.save("${npy}", emb)
print(f"Wrote ${npy} shape={emb.shape} key={key}")
PY

  _save_he_overview_if_needed "$sid"
  if [[ "${DELETE_SVS:-0}" == "1" ]]; then
    rm -f "$svs"
  fi
  if [[ "${DELETE_TILES:-1}" == "1" ]]; then
    rm -rf "$tiles_dir"
    rm -f "$tiles_h5ad"
  fi
  echo "[HE] done: ${sid}"
}

_submit_he_array() {
  if [[ ! -f "$IDS" ]]; then
    echo "Missing ${IDS}; run: $0 select" >&2
    exit 1
  fi
  local n max_idx array_spec job_script he_job
  n="$(grep -c . "$IDS" || true)"
  if [[ "$n" -lt 1 ]]; then
    echo "No sample ids in $IDS" >&2
    exit 1
  fi
  max_idx=$((n - 1))
  array_spec="0-${max_idx}"
  if [[ -n "${ARRAY_THROTTLE}" ]]; then
    array_spec="${array_spec}%${ARRAY_THROTTLE}"
  fi
  echo "Submitting HE GPU array ${array_spec} ($n slides) partition=${GPU_PARTITION}"
  echo "  fetch=${CONDA_ENV_NAME} tile=${CONDA_ENV_TILE_NAME} univ2=${CONDA_ENV_HE_NAME}"
  echo "  logs under ${LOG_ROOT}/<jobid>/"

  job_script="${LOG_ROOT}/_worker_he.sh"
  {
    cat <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=tcga_brca_he
#SBATCH --partition=${GPU_PARTITION}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=${HE_CPUS}
#SBATCH --mem=${HE_MEM}
#SBATCH --time=${HE_TIME}
#SBATCH --array=${array_spec}
set -euo pipefail
EOF
    _conda_bootstrap_snippet
    # Embed HE helpers + constants into the worker.
    declare -f _save_he_overview_if_needed _run_he_one_sample
    cat <<EOF
cd "${REPO_ROOT}"
SCRIPT_DIR="${SCRIPT_DIR}"
OUT_DIR="${OUT_DIR}"
SEL="${SEL}"
IDS="${IDS}"
PY_SCRIPT="${PY_SCRIPT}"
TILE_SCRIPT="${TILE_SCRIPT}"
UNIV2_SCRIPT="${UNIV2_SCRIPT}"
MODEL_DIR_RESOLVED="${MODEL_DIR_RESOLVED}"
CONDA_ENV_NAME="${CONDA_ENV_NAME}"
CONDA_ENV_TILE_NAME="${CONDA_ENV_TILE_NAME}"
CONDA_ENV_HE_NAME="${CONDA_ENV_HE_NAME}"
CROP_SIZE="${CROP_SIZE}"
DELETE_TILES="${DELETE_TILES:-1}"
DELETE_SVS="${DELETE_SVS:-0}"
HE_BATCH_SIZE="${HE_BATCH_SIZE:-32}"
FORCE_HE="${FORCE_HE:-0}"
FORCE_RETILE="${FORCE_RETILE:-0}"
IDX="\${SLURM_ARRAY_TASK_ID}"
SID="\$(sed -n "\$((IDX + 1))p" "\${IDS}")"
if [[ -z "\${SID}" ]]; then
  echo "Empty sample id at index \${IDX}" >&2
  exit 1
fi
echo "[HE] index=\${IDX} sample=\${SID} crop=\${CROP_SIZE}"
_run_he_one_sample "\${SID}"
EOF
  } > "$job_script"
  chmod +x "$job_script"

  # Hold → mkdir logs/<jobid> → release so %A path exists before tasks open logs.
  he_job="$(sbatch --parsable --hold \
    --output="${LOG_ROOT}/%A/he_%a.out" \
    --error="${LOG_ROOT}/%A/he_%a.err" \
    "$job_script")"
  he_job="${he_job%%;*}"
  mkdir -p "${LOG_ROOT}/${he_job}"
  scontrol release "${he_job}"
  echo "Submitted batch job ${he_job}"
  echo "HE logs: ${LOG_ROOT}/${he_job}/"
  printf '%s\n' "$he_job"
}

_submit_mudata() {
  local dep="${1:-}"
  local log_dir job_script mudata_job
  if [[ -n "$dep" ]]; then
    log_dir="${LOG_ROOT}/${dep}"
  else
    log_dir=""
  fi
  job_script="${LOG_ROOT}/_worker_mudata.sh"
  {
    cat <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=tcga_brca_mudata
#SBATCH --partition=${CPU_PARTITION}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${CPU_MEM}
#SBATCH --time=${CPU_TIME}
set -euo pipefail
EOF
    _conda_bootstrap_snippet
    cat <<EOF
_deepspotm_conda_activate ${CONDA_ENV_NAME}
cd "${REPO_ROOT}"
echo "[mudata] assembling MuData (all genes; no PCC)"
python "${PY_SCRIPT}" \\
  --out-dir "${OUT_DIR}" \\
  --stage mudata \\
  --he-encoder univ2 \\
  --univ2-model-dir "${MODEL_DIR_RESOLVED}"
EOF
  } > "$job_script"
  chmod +x "$job_script"

  if [[ -n "$dep" ]]; then
    mkdir -p "$log_dir"
    mudata_job="$(sbatch --parsable --hold --dependency=afterok:"${dep}" \
      --output="${log_dir}/mudata_%j.out" \
      --error="${log_dir}/mudata_%j.err" \
      "$job_script")"
    mudata_job="${mudata_job%%;*}"
    scontrol release "${mudata_job}"
    echo "Submitted batch job ${mudata_job} (afterok:${dep})"
    echo "Mudata logs: ${log_dir}/"
  else
    mudata_job="$(sbatch --parsable --hold \
      --output="${LOG_ROOT}/%j/mudata.out" \
      --error="${LOG_ROOT}/%j/mudata.err" \
      "$job_script")"
    mudata_job="${mudata_job%%;*}"
    mkdir -p "${LOG_ROOT}/${mudata_job}"
    scontrol release "${mudata_job}"
    echo "Submitted batch job ${mudata_job}"
    echo "Mudata logs: ${LOG_ROOT}/${mudata_job}/"
  fi
}

_run_he_test_one() {
  local sid="${1:-}"
  if [[ ! -f "$IDS" ]]; then
    echo "Missing ${IDS}; run: $0 select" >&2
    exit 1
  fi
  if [[ -z "$sid" ]]; then
    sid="$(sed -n '1p' "$IDS")"
  fi
  if [[ -z "$sid" ]]; then
    echo "No sample id (pass SAMPLE_ID or populate ${IDS})" >&2
    exit 1
  fi
  if ! grep -qxF "$sid" "$IDS"; then
    echo "Sample ${sid} not in ${IDS}" >&2
    exit 1
  fi
  export DELETE_TILES="${DELETE_TILES:-0}"
  export FORCE_HE="${FORCE_HE:-1}"
  export FORCE_RETILE="${FORCE_RETILE:-0}"
  export HE_BATCH_SIZE="${HE_BATCH_SIZE:-32}"
  echo "[test] HE pipeline locally for ${sid}"
  echo "  DELETE_TILES=${DELETE_TILES} FORCE_HE=${FORCE_HE} FORCE_RETILE=${FORCE_RETILE} CROP_SIZE=${CROP_SIZE}"
  echo "  fetch=${CONDA_ENV_NAME} tile=${CONDA_ENV_TILE_NAME} univ2=${CONDA_ENV_HE_NAME}"
  echo "  reuses tiles/{id}.h5ad if present; needs GPU for UNI2"
  _run_he_one_sample "$sid"
}

case "$STAGE" in
  select)
    _run_select_local
    ;;
  he)
    _submit_he_array
    ;;
  mudata)
    _submit_mudata
    ;;
  test)
    TEST_SID=""
    if ((${#USER_EXTRA[@]})) && [[ "${USER_EXTRA[0]}" != -* ]]; then
      TEST_SID="${USER_EXTRA[0]}"
    fi
    _run_he_test_one "$TEST_SID"
    ;;
  all)
    _run_select_local
    HE_JOB="$(_submit_he_array | tee /dev/stderr | awk '/^Submitted batch job / {print $4}' | tail -1)"
    if [[ -n "${HE_JOB}" ]]; then
      echo "HE array job: ${HE_JOB}; submitting mudata afterok"
      _submit_mudata "$HE_JOB"
    else
      echo "Could not parse HE job id; submit mudata manually after HE finishes:" >&2
      echo "  $0 mudata" >&2
    fi
    ;;
esac

echo "Log root: ${LOG_ROOT}/<jobid>/"
echo "Sample list: ${IDS}"
