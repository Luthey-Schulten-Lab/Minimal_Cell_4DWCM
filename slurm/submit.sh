#!/usr/bin/env bash
###############################################################################
# Convenience wrapper around slurm/4dwcm.sbatch.
#
#   ./slurm/submit.sh -o run001 -t 6300
#   ./slurm/submit.sh -o run002 -t 6300 -s 7 --time 24:00:00
#   ./slurm/submit.sh -o run001 -t 12600 --restart
#
# Every run reserves 2 GPUs (RDME solver + DNA/LAMMPS step).
###############################################################################
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

OUTDIR=""; SIM_TIME=""; SEED=1; MEMBRANE=1
IMAGE="4dwcm:ampere"; WALLTIME=""; DATA_ROOT="${REPO}/runs"
RESTART=0; MOUNT_REPO=0; DRYRUN=0

usage() {
  cat <<'EOF'
Usage: submit.sh -o OUTDIR -t SIM_TIME [options]

Required:
  -o, --outdir NAME     run name; output -> <data-root>/NAME/
  -t, --time-sim SECS   biological simulation time in seconds

Options:
  -s, --seed N          DNA RNG seed modifier      (default 1)
  -m, --membrane N      membrane flag              (default 1)
  -i, --image TAG       docker image               (default 4dwcm:ampere)
      --data-root DIR   output root                (default <repo>/runs)
      --time HH:MM:SS   SLURM wall limit; the sim stops ~20 min early to
                        checkpoint cleanly. Omit for no limit.
      --restart         resume via Restart_Whole_Cell_Minimal_Cell.py
      --mount-repo      run the working tree instead of the code baked
                        into the image (for iterating without rebuilding)
      --dry-run         print the sbatch command without submitting
  -h, --help            this message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--outdir)     OUTDIR="$2"; shift 2 ;;
    -t|--time-sim)   SIM_TIME="$2"; shift 2 ;;
    -s|--seed)       SEED="$2"; shift 2 ;;
    -m|--membrane)   MEMBRANE="$2"; shift 2 ;;
    -i|--image)      IMAGE="$2"; shift 2 ;;
    --data-root)     DATA_ROOT="$2"; shift 2 ;;
    --time)          WALLTIME="$2"; shift 2 ;;
    --restart)       RESTART=1; shift ;;
    --mount-repo)    MOUNT_REPO=1; shift ;;
    --dry-run)       DRYRUN=1; shift ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${OUTDIR}"   ]] && { echo "ERROR: -o/--outdir is required" >&2; usage; exit 1; }
[[ -z "${SIM_TIME}" ]] && { echo "ERROR: -t/--time-sim is required" >&2; usage; exit 1; }

# Fail early rather than after the job sits in the queue.
if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
  echo "ERROR: docker image '${IMAGE}' not found. Build it with ./docker/build_ampere.sh" >&2
  exit 1
fi

if [[ "${RESTART}" == "1" && ! -f "${DATA_ROOT}/${OUTDIR}/sim_properties.pkl" ]]; then
  echo "ERROR: --restart needs ${DATA_ROOT}/${OUTDIR}/sim_properties.pkl" >&2
  echo "       That run hasn't produced a checkpoint yet." >&2
  exit 1
fi

mkdir -p "${REPO}/slurm/logs" "${DATA_ROOT}"

SB=(sbatch --job-name "4dwcm-${OUTDIR}")
[[ -n "${WALLTIME}" ]] && SB+=(--time "${WALLTIME}")
SB+=(
  --export="ALL,REPO=${REPO},IMAGE=${IMAGE},OUTDIR=${OUTDIR},SIM_TIME=${SIM_TIME},SEED=${SEED},MEMBRANE=${MEMBRANE},DATA_ROOT=${DATA_ROOT},RESTART=${RESTART},MOUNT_REPO=${MOUNT_REPO}"
  "${REPO}/slurm/4dwcm.sbatch"
)

if [[ "${DRYRUN}" == "1" ]]; then
  printf '%q ' "${SB[@]}"; echo
  exit 0
fi

"${SB[@]}"
echo "output will appear in ${DATA_ROOT}/${OUTDIR}/"
echo "logs: ${REPO}/slurm/logs/"
