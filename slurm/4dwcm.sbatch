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
# No --partition is set, so the job lands on the cluster default. Override with
# `submit.sh -p NAME` or `sbatch -p NAME`.
#
# Submit via slurm/submit.sh (handles arg checking), or directly:
#   OUTDIR=run001 SIM_TIME=6300 sbatch slurm/4dwcm.sbatch
###############################################################################
#SBATCH --job-name=4dwcm
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

# ----------------------------------------------------------- CPU allocation
# `docker run` talks to the daemon, so the container is NOT placed in this job's
# cgroup and --cpus-per-task alone does not confine it. Pass the job's cpuset to
# Docker explicitly, otherwise concurrent runs land on the same physical cores.
#
# Within that set the RDME host thread and the btree_chromo BD process run at
# the same time and contend, so split it: LM takes the first half (applied below
# via taskset), btree_chromo the second (SpatialDnaDynamics reads DNA_CPU_CORES).
# Set either variable by hand to override. If the cpuset cannot be read, or is
# too small to divide, nothing is pinned.
#
# Split by *physical* core, not by logical cpu id. SLURM hands out hyperthread
# siblings together, so a sorted allocation looks like "8-15,120-127" where the
# second range is the sibling set of the first. Halving that list would put LM
# and btree_chromo on the same physical cores and defeat the isolation.
CPUSET=""
if command -v taskset >/dev/null 2>&1; then
  CPUSET="$(taskset -cp $$ 2>/dev/null | sed 's/.*:[[:space:]]*//')"
fi
if [[ -n "${CPUSET}" && ( -z "${LM_CPU_CORES:-}" || -z "${DNA_CPU_CORES:-}" ) ]]; then
  # "8-15,120-127" -> one cpu id per line
  mapfile -t CPUS < <(awk -F, '{
    for (i = 1; i <= NF; i++) {
      if (split($i, r, "-") == 2) { for (c = r[1]; c <= r[2]; c++) print c }
      else print $i
    }
  }' <<<"${CPUSET}")

  # One line per physical core: its siblings that are inside our cpuset.
  declare -A IN_SET=() && for c in "${CPUS[@]}"; do IN_SET["${c}"]=1; done
  CORES=() && declare -A SEEN=()
  for c in "${CPUS[@]}"; do
    [[ -n "${SEEN[${c}]:-}" ]] && continue
    SIB_FILE="/sys/devices/system/cpu/cpu${c}/topology/thread_siblings_list"
    SIBS=("${c}")
    if [[ -r "${SIB_FILE}" ]]; then
      mapfile -t ALL_SIBS < <(awk -F, '{
        for (i = 1; i <= NF; i++) {
          if (split($i, r, "-") == 2) { for (j = r[1]; j <= r[2]; j++) print j }
          else print $i
        }
      }' "${SIB_FILE}")
      SIBS=()
      for s in "${ALL_SIBS[@]}"; do [[ -n "${IN_SET[${s}]:-}" ]] && SIBS+=("${s}"); done
    fi
    for s in "${SIBS[@]}"; do SEEN["${s}"]=1; done
    CORES+=("$(IFS=,; echo "${SIBS[*]}")")
  done

  if (( ${#CORES[@]} >= 2 )); then
    HALF=$(( ${#CORES[@]} / 2 ))
    LM_CPU_CORES="${LM_CPU_CORES:-$(IFS=,; echo "${CORES[*]:0:HALF}")}"
    DNA_CPU_CORES="${DNA_CPU_CORES:-$(IFS=,; echo "${CORES[*]:HALF}")}"
  fi
fi

CPUSET_ARG=()
[[ -n "${CPUSET}" ]] && CPUSET_ARG=(--cpuset-cpus "${CPUSET}")

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

CMD=(python "${SCRIPT}"
     -od "${OUTDIR}"
     -t "${SIM_TIME}"
     -cd 0
     -drs "${SEED}"
     -m "${MEMBRANE}"
     -dsd /Software/
     "${MH_ARG[@]}")
[[ -n "${LM_CPU_CORES:-}" ]] && CMD=(taskset -c "${LM_CPU_CORES}" "${CMD[@]}")

echo "=== 4DWCM run ==============================================="
echo "  job          ${SLURM_JOB_ID:-<no slurm>} on $(hostname)"
echo "  image        ${IMAGE}"
echo "  script       ${SCRIPT}"
echo "  output       ${DATA_ROOT}/${OUTDIR}/"
echo "  sim time     ${SIM_TIME} s"
echo "  seed         ${SEED}"
echo "  host GPUs    ${GPUS}  -> container 0 (RDME) + 1 (DNA)"
echo "  host CPUs    ${CPUSET:-<unconfined>}"
echo "  cpu split    LM=${LM_CPU_CORES:-<unpinned>}  btree=${DNA_CPU_CORES:-<unpinned>}"
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
  -e DNA_CPU_CORES="${DNA_CPU_CORES:-}" \
  -e USER="$(id -un)" \
  -e LOGNAME="$(id -un)" \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  "${CPUSET_ARG[@]}" \
  "${MOUNTS[@]}" \
  "${IMAGE}" \
  "${CMD[@]}"

echo "=== finished: ${DATA_ROOT}/${OUTDIR}/ ==="
