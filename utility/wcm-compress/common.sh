# Shared helpers for pack_run.sh / extract_run.sh (sourced).
#   WCM_IMAGE   docker image to run the Python tools in (e.g. 4dwcm:ampere, 4dwcm:blackwell); empty = run them natively in the
#               current environment (python3 with numpy, and the zstd CLI >= 1.4 on PATH)
#   CORES       optional cpuset to pin to (e.g. 0-27); JOBS = worker processes (default: number of usable cores)
#   WCM_GPU_LOCKDIR  optional directory of gpu<N>.lock files, flock'ed before a GPU is used (shared-node etiquette)
A=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")"); WCM_IMAGE=${WCM_IMAGE:-}; CORES=${CORES:-}
JOBS=${JOBS:-$(python3 -c "import os;print(len(os.sched_getaffinity(0)) if not '$CORES' else sum(len(range(int(a),int(b or a)+1)) for a,_,b in (p.partition('-') for p in '$CORES'.split(','))))")}
pin() { if [ -n "$CORES" ]; then taskset -c "$CORES" "$@"; else "$@"; fi; }
# run <paths to bind...> -- <shell command>: in $WCM_IMAGE (each path bound at the same path) or natively
run() { local M=(); while [ "$1" != "--" ]; do M+=(-v "$1:$1"); shift; done; shift
  if [ -z "$WCM_IMAGE" ]; then pin bash -c "$*"; return; fi
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp ${CORES:+--cpuset-cpus "$CORES"} -v "$A:$A:ro" "${M[@]}" \
    -e PYTHONDONTWRITEBYTECODE=1 --entrypoint bash "$WCM_IMAGE" \
    -c "[ -f /opt/conda/etc/profile.d/conda.sh ] && source /opt/conda/etc/profile.d/conda.sh && conda activate \${CONDA_ENV:-lm_2.5_dev} 2>/dev/null; $*"; }
