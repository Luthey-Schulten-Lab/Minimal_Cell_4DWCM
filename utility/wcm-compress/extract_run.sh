#!/bin/bash
# Archive → run directory.
#   usage: extract_run.sh <archive_dir> <out_dir> [--keep-sidecar] [--skip-restart-files]
#   env: see common.sh; GPUS (GPU ids to try for the restart files, default: all), MAX_GPUS (3), JOBS_PER_GPU (24)
# 1 unpack (py/wcmpack.py) — for a full-layout run this is the whole job.
# 2 lean run (wcm_sidecar/LEAN.json, written by the lean-output writer): regenerate what the run did not write, byte for byte —
#   DNA/chromosome.lammpstrj (py/mono2lammpstrj.py + the literal frames in the sidecar), counts_fluxes_temp/*.csv (py/cftemp_regen.py),
#   DNA/data.lammps_<step> (btree_chromo, from the build carried in the sidecar, replays each hook's own directive file; needs GPUs;
#   the first hook's file comes from the sidecar) — then modes/mtimes, and wcm_sidecar/ is removed.
set -u; source "$(dirname "$(readlink -f "$0")")/common.sh"
ARCH=$(realpath "$1"); OUT=$(realpath -m "$2"); shift 2; KEEP=""; SKIPR=0
while [ $# -gt 0 ]; do case "$1" in --keep-sidecar) KEEP=--keep-sidecar;; --skip-restart-files) SKIPR=1;; esac; shift; done
t0=$(date +%s); [ -e "$OUT" ] && { echo "[extract] $OUT exists"; exit 1; }; mkdir -p "$(dirname "$OUT")"
el() { echo "$(( $(date +%s) - t0 )) s"; }
echo "[extract] $ARCH → $OUT (${WCM_IMAGE:-native}, ${JOBS} jobs)"
run "$ARCH" "$(dirname "$OUT")" -- "python3 $A/py/wcmpack.py unpack $ARCH $OUT --jobs $JOBS" > /dev/null || { echo "[extract] unpack failed"; exit 2; }
echo "[extract] unpacked ($(el))"
SC=$OUT/wcm_sidecar
if [ ! -f "$SC/LEAN.json" ]; then echo "[extract] full-layout run: done ($(el))"; exit 0; fi
NAME=$(python3 -c "import json;print(json.load(open('$SC/LEAN.json'))['run_dir_name'])")
BUILD=$(python3 -c "import json;print(json.load(open('$SC/LEAN.json')).get('btree_build_dir','build_fused'))")
python3 "$A/py/finish_extract.py" pre "$OUT"
run "$OUT" -- "cd $A/py && python3 mono2lammpstrj.py $OUT -o $OUT/DNA/chromosome.lammpstrj --literal $SC/literal --jobs $JOBS" | tail -1
echo "[extract] trajectory ($(el))"
pin python3 "$A/py/cftemp_regen.py" "$OUT" "$OUT/counts_fluxes_temp" | tail -1
if [ $SKIPR = 0 ]; then
  cp -p "$SC"/literal/DNA/data.lammps_* "$OUT/DNA/"; FIRST=$(ls "$SC/literal/DNA/" | sed 's/data.lammps_//' | paste -sd,)
  ALL=$(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null | tr '\n' ' '); GL=(); n=0
  for g in ${GPUS:-$ALL}; do [ $n -ge ${MAX_GPUS:-3} ] && break
    if [ -n "${WCM_GPU_LOCKDIR:-}" ]; then exec {fd}>"$WCM_GPU_LOCKDIR/gpu$g.lock"; flock -n $fd || continue; fi
    GL+=($g); n=$((n+1)); done
  [ ${#GL[@]} = 0 ] && { echo "[extract] no GPU for the restart files (set GPUS, or --skip-restart-files)"; exit 3; }
  SCR=$OUT.replay_tmp; mkdir -p "$SCR"; i=0
  for g in "${GL[@]}"; do
    ARGS="$NAME /out part:$i/${#GL[@]} --jobs ${JOBS_PER_GPU:-24} --emit /rundir/DNA --skip $FIRST"
    if [ -n "$WCM_IMAGE" ]; then
      docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp --gpus device=$g ${CORES:+--cpuset-cpus "$CORES"} -e BTREE_BUILD_DIR=$BUILD \
        -v "$A:$A:ro" -v "$OUT:/rundir" -v "$SC/btree:/Software/btree_chromo:ro" -v "$SC/input_data:/inp:ro" -v "$SCR:/out" --entrypoint bash "$WCM_IMAGE" \
        -c "[ -f /opt/conda/etc/profile.d/conda.sh ] && source /opt/conda/etc/profile.d/conda.sh && conda activate \${CONDA_ENV:-lm_2.5_dev} 2>/dev/null; python3 $A/py/data_replay.py /rundir $ARGS" > "$SCR/part$i.txt" 2>&1 &
    else
      mkdir -p "$SCR/p$i"; ( cd "$SCR/p$i" && CUDA_VISIBLE_DEVICES=$g BTREE_DIR="$SC/btree" INPUT_DIR="$SC/input_data" BTREE_BUILD_DIR=$BUILD \
        pin python3 "$A/py/data_replay.py" "$OUT" $(echo "$ARGS" | sed "s#/out#$SCR/p$i#; s#/rundir/DNA#$OUT/DNA#") ) > "$SCR/part$i.txt" 2>&1 &
    fi
    i=$((i+1)); done; wait
  echo "[extract] restart files on GPUs ${GL[*]} ($(el)): $(grep -h '"steps"' "$SCR"/part*.txt | paste -sd' ')"
  for try in 1 2 3; do                                             # a replay can fail transiently under load: redo the missing ones
    MISS=$(python3 -c "import os,re;D='$OUT/DNA';print(','.join(m.group(1) for m in (re.match(r'dna_monomers_(\\d+)\\.bin$',f) for f in sorted(os.listdir(D))) if m and not os.path.exists(D+'/data.lammps_'+m.group(1))))")
    [ -z "$MISS" ] && break; echo "[extract] retry $try: $MISS"; g=${GL[0]}
    RARGS="$NAME /out $MISS --jobs 4 --emit /rundir/DNA"
    if [ -n "$WCM_IMAGE" ]; then
      docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp --gpus device=$g ${CORES:+--cpuset-cpus "$CORES"} -e BTREE_BUILD_DIR=$BUILD \
        -v "$A:$A:ro" -v "$OUT:/rundir" -v "$SC/btree:/Software/btree_chromo:ro" -v "$SC/input_data:/inp:ro" -v "$SCR:/out" --entrypoint bash "$WCM_IMAGE" \
        -c "[ -f /opt/conda/etc/profile.d/conda.sh ] && source /opt/conda/etc/profile.d/conda.sh && conda activate \${CONDA_ENV:-lm_2.5_dev} 2>/dev/null; python3 $A/py/data_replay.py /rundir $RARGS" >> "$SCR/retry.txt" 2>&1
    else
      mkdir -p "$SCR/r"; ( cd "$SCR/r" && CUDA_VISIBLE_DEVICES=$g BTREE_DIR="$SC/btree" INPUT_DIR="$SC/input_data" BTREE_BUILD_DIR=$BUILD \
        pin python3 "$A/py/data_replay.py" "$OUT" $(echo "$RARGS" | sed "s#/out#$SCR/r#; s#/rundir/DNA#$OUT/DNA#") ) >> "$SCR/retry.txt" 2>&1
    fi
  done
  rm -rf "$SCR"
fi
python3 "$A/py/finish_extract.py" post "$OUT" $KEEP $([ $SKIPR = 1 ] && echo --skip-restart-files); rc=$?
echo "[extract] done ($(el))"; exit $rc
