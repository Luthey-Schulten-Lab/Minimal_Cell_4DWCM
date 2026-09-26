#!/usr/bin/env python3
"""Lossless check: the unpacked tree must equal the manifest of the original run — every file's sha256, size, mode, mtime_ns;
the same directory set (and modes); the same symlinks. usage: gate.py <manifest.json> <unpacked_dir> [nproc] → exit 0 = PASS

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import sys, json, os, subprocess, tempfile
base = json.load(open(sys.argv[1])); root = sys.argv[2]; n = sys.argv[3] if len(sys.argv) > 3 else '24'
tmp = tempfile.mktemp(suffix='.json')
subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), 'manifest.py'), root, tmp, n], check=True, stdout=subprocess.DEVNULL)
got = json.load(open(tmp)); os.remove(tmp)
bad = []
B = {e['path']: e for e in base['files']}; G = {e['path']: e for e in got['files']}
for p in sorted(set(B) | set(G)):
    if p not in G: bad.append('missing ' + p); continue
    if p not in B: bad.append('extra ' + p); continue
    for k in ('sha256', 'size', 'mode', 'mtime_ns'):
        if B[p][k] != G[p][k]: bad.append('%s %s: %s != %s' % (k, p, G[p][k], B[p][k]))
bd = {d['path']: d['mode'] for d in base['dirs']}; gd = {d['path']: d['mode'] for d in got['dirs']}
if bd != gd: bad.append('dirs differ: %s' % sorted(set(bd.items()) ^ set(gd.items()))[:10])
if sorted(map(str, base['links'])) != sorted(map(str, got['links'])): bad.append('links differ')
print(json.dumps({'gate': 'PASS' if not bad else 'FAIL', 'files': len(B), 'bytes': base['bytes'], 'problems': bad[:30], 'n_problems': len(bad)}))
sys.exit(0 if not bad else 1)
