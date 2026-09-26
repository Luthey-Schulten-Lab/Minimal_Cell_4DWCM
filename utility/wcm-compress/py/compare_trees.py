#!/usr/bin/env python3
"""Compare an extracted run with a reference full run (e.g. a run made with WCM_LEAN_SHADOW=1, which writes both layouts). usage: compare_trees.py <reference> <extracted> [lean_view]
Every reference file (outside wcm_sidecar/) must exist with equal sha256, size and mode; files present in the lean view (written by the
lean run) must also keep their mtime. Extra files in the extracted tree are reported.

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, sys, json, hashlib
from concurrent.futures import ProcessPoolExecutor
ref, ext = sys.argv[1], sys.argv[2]; lean = sys.argv[3] if len(sys.argv) > 3 else None
def files(root):
    out = {}
    for d, ds, fs in os.walk(root):
        if os.path.relpath(d, root).split(os.sep)[0] == 'wcm_sidecar': continue
        for f in fs: out[os.path.relpath(os.path.join(d, f), root)] = os.path.join(d, f)
    return out
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''): h.update(b)
    return h.hexdigest()
R, X = files(ref), files(ext); keep = set(files(lean)) if lean else set()
common = sorted(set(R) & set(X)); bad = []
with ProcessPoolExecutor(28) as ex:
    hr = dict(zip(common, ex.map(sha, [R[p] for p in common], chunksize=16))); hx = dict(zip(common, ex.map(sha, [X[p] for p in common], chunksize=16)))
for p in common:
    a, b = os.stat(R[p]), os.stat(X[p])
    if hr[p] != hx[p]: bad.append('content ' + p)
    elif a.st_size != b.st_size or (a.st_mode & 0o7777) != (b.st_mode & 0o7777): bad.append('size/mode ' + p)
    elif p in keep and a.st_mtime_ns != b.st_mtime_ns: bad.append('mtime ' + p)
missing = sorted(set(R) - set(X)); extra = sorted(set(X) - set(R))
print(json.dumps({'reference_files': len(R), 'extracted_files': len(X), 'identical_content': len(common) - sum(b.startswith('content') for b in bad),
                  'regenerated_checked': len([p for p in common if p not in keep]), 'problems': bad[:20], 'n_problems': len(bad),
                  'missing': missing[:10], 'n_missing': len(missing), 'extra': extra[:10], 'n_extra': len(extra),
                  'result': 'IDENTICAL' if not bad and not missing and not extra else 'DIFFERENT'}))
