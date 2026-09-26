#!/bin/bash
#SBATCH --job-name=wcm-compress
#SBATCH --output=wcm-compress-%j.out
#SBATCH --time=0-03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
# Compress a finished run as a CPU-only job (submit it from the end of a run job so the GPUs are released right away).
#   RUNDIR=<run dir> [WCM_COMPRESS=<repo>/utility/wcm-compress] [WCM_IMAGE=...] sbatch utility/wcm-compress/pack_job.sbatch
# On a passing bit-for-bit check the run directory is replaced by <run dir>.wcmpack; otherwise it is left untouched.
set -u
WCM_COMPRESS=${WCM_COMPRESS:-${SLURM_SUBMIT_DIR:-.}/utility/wcm-compress}; RUNDIR=${RUNDIR:?set RUNDIR}
export CORES=${CORES:-$(taskset -cp $$ | sed 's/.*:[[:space:]]*//')}
echo "[pack-job] ${SLURM_JOB_ID:-none} $(date -Is) $RUNDIR"
"$WCM_COMPRESS/pack_run.sh" "$RUNDIR"; rc=$?
echo "[pack-job] rc=$rc $(date -Is)"; exit $rc
