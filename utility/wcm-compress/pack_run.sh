#!/bin/bash
# Compress one 4DWCM run directory losslessly, prove it, and replace the run with the archive.
#   usage: pack_run.sh <run_dir> [archive_dir]          (env: see common.sh; default archive: <run_dir>.wcmpack)
# 1 byte manifest of the run (sha256, size, mode, mtime of every file) → 2 pack (py/wcmpack.py) → 3 unpack into a scratch dir next to
# the archive → 4 check: every file, directory and link of the unpacked copy identical to the run, bit for bit
# → PASS: <archive>/VERIFIED.json (+ a copy of the decoder) and the run directory is deleted
# → FAIL (or any error): the archive is deleted and the run directory is left exactly as it was.
set -u; source "$(dirname "$(readlink -f "$0")")/common.sh"
RUN=$(realpath "$1"); ARCH=$(realpath -m "${2:-${RUN%/}.wcmpack}")
case "$ARCH/" in "$RUN"/*) echo "[pack] the archive must not be inside the run directory"; exit 1;; esac
fail() { echo "[pack] $1 — archive removed, run directory kept"; rm -rf "$ARCH" "$SCR"; exit $2; }
SCR=${ARCH}.verify_tmp; t0=$(date +%s); mkdir -p "$ARCH"; rm -rf "$SCR"
echo "[pack] $RUN → $ARCH (${WCM_IMAGE:-native}, ${JOBS} jobs${CORES:+, cores $CORES})"
pin python3 "$A/py/manifest.py" "$RUN" "$ARCH/source_manifest.json" "$JOBS" || fail "manifest failed" 2
run "$RUN" "$(dirname "$ARCH")" -- "python3 $A/py/wcmpack.py pack $RUN $ARCH --jobs $JOBS" > "$ARCH/pack.json" || fail "pack failed" 3
run "$(dirname "$ARCH")" -- "python3 $A/py/wcmpack.py unpack $ARCH $SCR --jobs $JOBS" > "$ARCH/unpack.json" 2>/dev/null || fail "unpack failed" 4
pin python3 "$A/py/gate.py" "$ARCH/source_manifest.json" "$SCR" "$JOBS" > "$ARCH/gate.json"; G=$?; rm -rf "$SCR"
cp "$A/py/wcmpack.py" "$ARCH/wcmpack.py"
python3 - "$ARCH" "$RUN" "$G" "$t0" <<'PY'
import json, sys, time, os
a, run, g, t0 = sys.argv[1:]; lj = lambda p: json.loads([l for l in open(p) if l.startswith('{')][-1])
p = lj(a + '/pack.json'); gt = lj(a + '/gate.json'); lean = os.path.exists(os.path.join(run, 'wcm_sidecar', 'LEAN.json'))
v = {'run': run, 'lean_run': lean, 'gate': gt['gate'], 'files': gt['files'], 'raw_bytes': p['raw_bytes'], 'packed_bytes': p['packed_bytes'],
     'ratio': round(p['raw_bytes'] / p['packed_bytes'], 3), 'wall_s': int(time.time() - float(t0)),
     'verified_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'extract': 'utility/wcm-compress/extract_run.sh %s <out_dir>' % a}
json.dump(v, open(a + '/VERIFIED.json', 'w'), indent=1)
print('[pack] %s  %.2f GB -> %.2f GB (%.2fx)  check %s  %d s' % (os.path.basename(run), v['raw_bytes'] / 1e9, v['packed_bytes'] / 1e9, v['ratio'], v['gate'], v['wall_s']))
PY
[ $G -ne 0 ] && { cp "$ARCH/gate.json" "${ARCH}.FAILED_check.json" 2>/dev/null; fail "CHECK FAILED (details: ${ARCH}.FAILED_check.json)" 5; }
# the run is still byte-identical to what was checked? (nothing may have written into it meanwhile)
pin python3 "$A/py/gate.py" "$ARCH/source_manifest.json" "$RUN" "$JOBS" > /dev/null || fail "the run directory changed while it was being packed" 6
rm -rf "$RUN" && echo "[pack] verified bit for bit; run directory replaced by $ARCH"; exit 0
