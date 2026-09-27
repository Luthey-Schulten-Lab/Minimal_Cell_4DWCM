#!/usr/bin/env bash
###############################################################################
# Submit one 4DWCM run from a private snapshot of this working tree, so nothing a run writes lands in the repo.
#
#   slurm/submit_run.sh -b perturb_20260927 -n q_ko_rnap -t 30 -g 2,3 [-p perturbations/tests/q_ko_rnap.yaml] [-s 1]
#
# Layout (ROOT defaults to /raid/racda/4dwcm_tests):
#   $ROOT/<batch>/<name>/               the run's Data directory (counts_and_fluxes.csv, perturbation_applied.json, ...)
#   $ROOT/<batch>/_src/<name>/4dwcm/    the code it ran: a copy of the tracked sources + modelspec spec + perturbation files
#   $ROOT/<batch>/_src/<name>/PROVENANCE.txt   commit, branch, working-tree diff, LM / btree builds
#   $ROOT/<batch>/logs/<name>-<jobid>.{out,err}
#
# odecell writes and compiles its generated ODE code (cythonCompiledFunctions.pyx, setup_tmp.py, build/, pyxbld/) in the
# working directory, with the rate constants compiled in; one snapshot per run keeps runs from sharing or overwriting it.
# The wcm-speed Python needs the ERA LM and btree builds; LM_BUILD / BTREE_REPO / BTREE_BUILD_DIR default to the ones
# the Sep 26 lean runs used.
###############################################################################
set -euo pipefail
ROOT="${ROOT:-/raid/racda/4dwcm_tests}"
SNAP_DEFAULT=/raid/racda/4dwcm_runs/_src/Sep26_leanB_rng48
LM_BUILD="${LM_BUILD:-${SNAP_DEFAULT}/lm_build}"
BTREE_REPO="${BTREE_REPO:-${SNAP_DEFAULT}/btree}"
BTREE_BUILD_DIR="${BTREE_BUILD_DIR:-build_fused}"
IMAGE="${IMAGE:-4dwcm:blackwell}"
BATCH="" NAME="" SIM_TIME="" GPUS="" PERTURB="" SEED=1
while getopts "b:n:t:g:p:s:" o; do
  case "$o" in b) BATCH=$OPTARG;; n) NAME=$OPTARG;; t) SIM_TIME=$OPTARG;; g) GPUS=$OPTARG;; p) PERTURB=$OPTARG;; s) SEED=$OPTARG;;
    *) sed -n '3,6p' "$0"; exit 2;; esac
done
[[ -n "$BATCH" && -n "$NAME" && -n "$SIM_TIME" && -n "$GPUS" ]] || { sed -n '3,6p' "$0"; exit 2; }
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
B="$ROOT/$BATCH"; SRC="$B/_src/$NAME"
[[ -e "$B/$NAME" ]] && { echo "ERROR: $B/$NAME exists; pick another name (never overwrite a run)." >&2; exit 1; }
[[ -e "$SRC" ]] && { echo "ERROR: $SRC exists." >&2; exit 1; }
[[ -z "$PERTURB" || -f "$REPO/$PERTURB" ]] || { echo "ERROR: $PERTURB not found in $REPO." >&2; exit 1; }
[[ -f "$REPO/modelspec/out/model_spec.json" ]] || echo "note: no modelspec/out/model_spec.json; a perturbed run rebuilds the spec at start"
mkdir -p "$SRC/4dwcm" "$B/logs"

# tracked sources as they are in the working tree (uncommitted edits included), plus the spec and perturbation files
( cd "$REPO" && git ls-files -z --cached --others --exclude-standard -- \
    ':!Lattice_Microbes' ':!btree_chromo_gpu' ':!lammps' ':!odecell' ':!sc_chain_generation' ':!analysis' ':!Data' \
    ':!*.zip' ':!*.swp' ':!4DWCM-GUI.html' ':!slurm/logs' ':!build' ':!cythonCompiledFunctions.c' ':!*.so' ) | \
  ( cd "$REPO" && xargs -0 cp --parents -t "$SRC/4dwcm" )
if [[ -f "$REPO/modelspec/out/model_spec.json" ]]; then    # gitignored build output: the spec a perturbed run validates against
  mkdir -p "$SRC/4dwcm/modelspec/out"; cp "$REPO/modelspec/out/model_spec.json" "$SRC/4dwcm/modelspec/out/"
fi
{
  echo "4DWCM run $NAME (batch $BATCH), snapshot $(date -Is)"
  echo "code      : $(git -C "$REPO" rev-parse --abbrev-ref HEAD) $(git -C "$REPO" rev-parse --short HEAD) from $REPO"
  echo "perturb   : ${PERTURB:-none}"
  echo "LM build  : $LM_BUILD"
  echo "btree     : $BTREE_REPO ($BTREE_BUILD_DIR)"
  echo "image     : $IMAGE   GPUs $GPUS   seed $SEED   bio time $SIM_TIME s"
  echo; echo "working-tree changes vs HEAD:"; git -C "$REPO" status --short -- . ':!Lattice_Microbes' ':!btree_chromo_gpu' ':!lammps' ':!odecell' ':!sc_chain_generation'
} > "$SRC/PROVENANCE.txt"
git -C "$REPO" diff > "$SRC/working_tree.diff"

export REPO="$SRC/4dwcm" MOUNT_REPO=1 DATA_ROOT="$B" OUTDIR="$NAME" SIM_TIME SEED GPUS PERTURB IMAGE LM_BUILD BTREE_REPO BTREE_BUILD_DIR
sbatch --job-name="wcm_$NAME" --output="$B/logs/$NAME-%j.out" --error="$B/logs/$NAME-%j.err" --export=ALL "$REPO/slurm/4dwcm.sbatch"
