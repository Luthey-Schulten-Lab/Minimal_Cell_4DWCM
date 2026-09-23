#!/usr/bin/env bash
###############################################################################
# Submit one 4DWCM replicate from its own worktree.
#
#   ./slurm/submit_modularize.sh -s 51 -g 0,2
#   ./slurm/submit_modularize.sh -s 52 -g 3,4 --dry-run
#
# submit.sh derives REPO from its own location, so it always mounts the
# checkout it lives in. Production replicates each run from a separate worktree
# -- the container writes generated Cython sources into the tree it is given,
# so two concurrent runs cannot share one -- which means REPO has to be chosen
# per run instead. That is the whole reason this wrapper exists.
#
# The worktree for a seed is created from --ref if it is not already there, so
# re-running this for a finished seed reuses the tree rather than rebuilding it.
###############################################################################
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

SEED=""; GPUS=""; SIM_TIME=7200; OUTDIR=""; MEMBRANE=1; DRYRUN=0
IMAGE="${IMAGE:-4dwcm:blackwell}"
REF="${REF:-origin/modularize}"

# Worktrees land beside the checkout and output beside them; both are usually
# on a scratch filesystem rather than $HOME, hence the env overrides.
WORKTREE_ROOT="${WORKTREE_ROOT:-$(dirname "${REPO_ROOT}")}"
DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/runs}"

# Bind a btree_chromo checkout over the one in the image. Only needed while the
# image predates a btree fix you want; an image built from the current
# Dockerfile already carries the pinned BTREE_REF, so leave this empty then.
BTREE_REPO="${BTREE_REPO:-}"
BTREE_BUILD_DIR="${BTREE_BUILD_DIR:-build_fused}"

usage() {
  cat <<'EOF'
Usage: submit_modularize.sh -s SEED -g GPU,GPU [options]

Required:
  -s, --seed N          DNA RNG seed; also names the worktree and output dir
  -g, --gpus A,B        two host GPU indices: RDME card, then DNA card

Options:
  -t, --time-sim SECS   biological seconds            (default 7200)
  -o, --outdir NAME     output dir name               (default <MonDD>_modularize_rng<seed>)
  -m, --membrane N      membrane flag                 (default 1)
  -i, --image TAG       docker image, or $IMAGE       (default 4dwcm:blackwell)
      --ref REF         git ref for a new worktree    (default origin/modularize)
      --data-root DIR   output root, or $DATA_ROOT
      --worktree-root DIR  where worktrees live, or $WORKTREE_ROOT
      --btree DIR       btree_chromo checkout to bind over the image's copy
      --btree-build DIR build subdir within it        (default build_fused)
      --dry-run         print what would be submitted
  -h, --help            this message

Pick GPUs that are actually idle. The scheduler only tracks work it started, so
it will hand out cards that are busy with anything else on the box; 4dwcm.sbatch
refuses to start on those rather than running slowly next to them.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--seed)          SEED="$2"; shift 2 ;;
    -g|--gpus)          GPUS="$2"; shift 2 ;;
    -t|--time-sim)      SIM_TIME="$2"; shift 2 ;;
    -o|--outdir)        OUTDIR="$2"; shift 2 ;;
    -m|--membrane)      MEMBRANE="$2"; shift 2 ;;
    -i|--image)         IMAGE="$2"; shift 2 ;;
    --ref)              REF="$2"; shift 2 ;;
    --data-root)        DATA_ROOT="$2"; shift 2 ;;
    --worktree-root)    WORKTREE_ROOT="$2"; shift 2 ;;
    --btree)            BTREE_REPO="$2"; shift 2 ;;
    --btree-build)      BTREE_BUILD_DIR="$2"; shift 2 ;;
    --dry-run)          DRYRUN=1; shift ;;
    -h|--help)          usage; exit 0 ;;
    *) echo "unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "${SEED}" ]] && { echo "ERROR: -s/--seed is required" >&2; usage; exit 1; }
[[ -z "${GPUS}" ]] && { echo "ERROR: -g/--gpus is required, e.g. -g 0,2" >&2; usage; exit 1; }

if [[ ! "${SEED}" =~ ^[0-9]+$ ]]; then
  echo "ERROR: seed must be a number, got '${SEED}'" >&2; exit 1
fi
if [[ ! "${GPUS}" =~ ^[0-9]+,[0-9]+$ ]]; then
  echo "ERROR: --gpus needs exactly two indices, e.g. 0,2 (got '${GPUS}')" >&2; exit 1
fi

OUTDIR="${OUTDIR:-$(date +%b%d)_modularize_rng${SEED}}"
WORKTREE="${WORKTREE_ROOT}/modularize_run_rng${SEED}"
SBATCH="${REPO_ROOT}/slurm/4dwcm.sbatch"

[[ -x "${SBATCH}" ]] || { echo "ERROR: ${SBATCH} missing or not executable" >&2; exit 1; }

if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
  echo "ERROR: docker image '${IMAGE}' not found. Build one with docker/build_blackwell.sh" >&2
  exit 1
fi

if [[ -n "${BTREE_REPO}" && ! -x "${BTREE_REPO}/${BTREE_BUILD_DIR}/apps/btree_chromo" ]]; then
  echo "ERROR: --btree given but ${BTREE_REPO}/${BTREE_BUILD_DIR}/apps/btree_chromo is not executable" >&2
  exit 1
fi

if [[ ! -d "${WORKTREE}" ]]; then
  if [[ "${DRYRUN}" == "1" ]]; then
    echo "would create worktree ${WORKTREE} at ${REF}"
  else
    git -C "${REPO_ROOT}" worktree add --detach "${WORKTREE}" "${REF}"
  fi
fi

# 4dwcm.sbatch writes its log to a path relative to the submission directory,
# so submit from the worktree with that directory already in place. Otherwise
# SLURM silently discards the job output.
[[ "${DRYRUN}" == "1" ]] || mkdir -p "${WORKTREE}/slurm/logs" "${DATA_ROOT}"

if [[ "${DRYRUN}" == "1" ]]; then
  cat <<EOF
would submit:
  job        4dwcm-rng${SEED}
  worktree   ${WORKTREE}  (ref ${REF})
  image      ${IMAGE}
  output     ${DATA_ROOT}/${OUTDIR}/
  sim time   ${SIM_TIME} s
  gpus       ${GPUS}
  btree      ${BTREE_REPO:-<baked into image>}${BTREE_REPO:+ (${BTREE_BUILD_DIR})}
  logs       ${WORKTREE}/slurm/logs/
EOF
  exit 0
fi

# Export rather than listing these in --export: SLURM splits that list on
# commas, which would cut GPUS=0,2 in half and leave the job holding one card.
cd "${WORKTREE}"
export REPO="${WORKTREE}"
export IMAGE OUTDIR SIM_TIME SEED MEMBRANE DATA_ROOT GPUS
export BTREE_REPO BTREE_BUILD_DIR
export RESTART=0 MOUNT_REPO=1

sbatch --job-name="4dwcm-rng${SEED}" --export=ALL "${SBATCH}"
echo "output -> ${DATA_ROOT}/${OUTDIR}/"
echo "logs   -> ${WORKTREE}/slurm/logs/"
