#!/usr/bin/env bash
###############################################################################
# 4DWCM single simulation — SLURM + Docker, 2 GPUs per run.
#
# One run needs two GPUs:
#   RDME / Lattice Microbes solver  -> -cd 0          (container-local index)
#   btree_chromo / LAMMPS DNA step  -> DNA_GPU_ID=1   (container-local index)
#
# SLURM hands us physical GPUs in $CUDA_VISIBLE_DEVICES (e.g. "3,6"); those are
# passed to Docker, which renumbers them to 0 and 1 inside the container. So the
# indices above are always 0/1 regardless of which physical cards we were given.
#
# Submit via slurm/submit.sh (handles arg checking), or directly:
#   OUTDIR=run001 SIM_TIME=6300 sbatch slurm/4dwcm.sbatch
###############################################################################
#SBATCH --job-name=4dwcm
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/logs/4dwcm-%j.out
#SBATCH --error=slurm/logs/4dwcm-%j.err

set -euo pipefail

# ---------------------------------------------------------------- parameters
IMAGE="${IMAGE:-4dwcm:ampere}"
OUTDIR="${OUTDIR:?set OUTDIR to a run name, e.g. OUTDIR=run001}"
SIM_TIME="${SIM_TIME:?set SIM_TIME in biological seconds, e.g. SIM_TIME=6300}"
SEED="${SEED:-1}"
MEMBRANE="${MEMBRANE:-1}"
DATA_ROOT="${DATA_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/runs}"
MOUNT_REPO="${MOUNT_REPO:-0}"     # 1 = use the working tree instead of baked-in code
RESTART="${RESTART:-0}"           # 1 = resume an existing run via the Restart script

# SLURM copies this script into its spool dir, so BASH_SOURCE does NOT point at
# the repo. Prefer REPO exported by submit.sh, then the submit directory, and
# only fall back to BASH_SOURCE when run outside SLURM.
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"

# Output lands in $DATA_ROOT/$OUTDIR/ (bound to /src/4d/Data in the container).
mkdir -p "${DATA_ROOT}" "${REPO}/slurm/logs"

# ------------------------------------------------------------- GPU selection
# SLURM exports the allocated device indices; fall back to 0,1 outside SLURM.
GPUS="${CUDA_VISIBLE_DEVICES:-0,1}"
NGPU=$(awk -F, '{print NF}' <<<"${GPUS}")
if [[ "${NGPU}" -lt 2 ]]; then
  echo "ERROR: need 2 GPUs, SLURM gave ${NGPU} (CUDA_VISIBLE_DEVICES=${GPUS})." >&2
  echo "       Submit with --gres=gpu:2." >&2
  exit 1
fi

# ------------------------------------------------- wall-clock safety cut-off
# The sim stops cleanly at -mh so it can checkpoint before SLURM kills it.
# Derive it from the job's own time limit, minus a 20-minute margin.
MH_ARG=()
if [[ -n "${MAX_HOURS:-}" ]]; then
  MH_ARG=(-mh "${MAX_HOURS}")
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
  LIMIT=$(squeue -h -j "${SLURM_JOB_ID}" -o %l 2>/dev/null | tr -d ' ' || true)
  if [[ -n "${LIMIT}" && "${LIMIT}" != "UNLIMITED" && "${LIMIT}" != "NOT_SET" ]]; then
    SECS=$(awk -F'[-:]' '{
      if (NF==4)      print (($1*24+$2)*3600)+($3*60)+$4;
      else if (NF==3) print ($1*3600)+($2*60)+$3;
      else if (NF==2) print ($1*60)+$2;
      else            print $1*60;
    }' <<<"${LIMIT}")
    MARGIN=1200
    if [[ "${SECS}" -gt $((MARGIN + 600)) ]]; then
      MH_ARG=(-mh "$(awk -v s="${SECS}" -v m="${MARGIN}" 'BEGIN{printf "%.3f", (s-m)/3600}')")
    fi
  fi
fi

# ----------------------------------------------------------- container setup
SCRIPT="Whole_Cell_Minimal_Cell.py"
[[ "${RESTART}" == "1" ]] && SCRIPT="Restart_Whole_Cell_Minimal_Cell.py"

CIDFILE="$(mktemp -u "${TMPDIR:-/tmp}/4dwcm-cid-XXXXXX")"
cleanup() {
  if [[ -f "${CIDFILE}" ]]; then
    docker rm -f "$(cat "${CIDFILE}")" >/dev/null 2>&1 || true
    rm -f "${CIDFILE}"
  fi
}
# Without this, scancel kills the job step but leaves the container running.
trap cleanup EXIT TERM INT

# The code resolves every input as <workdir>/input_data/... and writes output to
# <workdir>/Data/<outdir>/. The image already holds code and input_data at
# /src/4d, so leave the workdir alone and bind only Data/ out to the host.
MOUNTS=(-v "${DATA_ROOT}:/src/4d/Data")
[[ "${MOUNT_REPO}" == "1" ]] && MOUNTS+=(-v "${REPO}:/src/4d")

echo "=== 4DWCM run ==============================================="
echo "  job          ${SLURM_JOB_ID:-<no slurm>} on $(hostname)"
echo "  image        ${IMAGE}"
echo "  script       ${SCRIPT}"
echo "  output       ${DATA_ROOT}/${OUTDIR}/"
echo "  sim time     ${SIM_TIME} s"
echo "  seed         ${SEED}"
echo "  host GPUs    ${GPUS}  -> container 0 (RDME) + 1 (DNA)"
echo "  wall cutoff  ${MH_ARG[*]:-none}"
echo "  repo mount   $([[ "${MOUNT_REPO}" == "1" ]] && echo "working tree" || echo "baked into image")"
echo "============================================================="

docker run --rm \
  --cidfile "${CIDFILE}" \
  --gpus "\"device=${GPUS}\"" \
  --user "$(id -u):$(id -g)" \
  --shm-size=8g \
  -e HOME=/tmp \
  -e DNA_GPU_ID=1 \
  -e USER="$(id -un)" \
  -e LOGNAME="$(id -un)" \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  "${MOUNTS[@]}" \
  "${IMAGE}" \
  python "${SCRIPT}" \
    -od "${OUTDIR}" \
    -t "${SIM_TIME}" \
    -cd 0 \
    -drs "${SEED}" \
    -m "${MEMBRANE}" \
    -dsd /Software/ \
    "${MH_ARG[@]}"

echo "=== finished: ${DATA_ROOT}/${OUTDIR}/ ==="
