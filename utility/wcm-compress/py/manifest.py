#!/usr/bin/env python3
"""Byte-level manifest of a run tree: every regular file (path, size, mode, mtime_ns, sha256), every dir, every symlink.
usage: manifest.py <root> <out.json> [nproc]   — the lossless gate compares two of these (tools/gate.py).

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, sys, json, hashlib, stat
from concurrent.futures import ProcessPoolExecutor
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda: f.read(1<<24), b''): h.update(b)
    return h.hexdigest()
def main(root, out, n=16):
    root=os.path.abspath(root); files=[]; dirs=[]; links=[]
    for d,ds,fs in os.walk(root):
        rd=os.path.relpath(d,root); st=os.lstat(d); dirs.append({'path':rd,'mode':stat.S_IMODE(st.st_mode)})
        for f in fs:
            p=os.path.join(d,f); st=os.lstat(p); rp=os.path.relpath(p,root)
            if stat.S_ISLNK(st.st_mode): links.append({'path':rp,'target':os.readlink(p)})
            elif stat.S_ISREG(st.st_mode): files.append({'path':rp,'size':st.st_size,'mode':stat.S_IMODE(st.st_mode),'mtime_ns':st.st_mtime_ns})
            else: raise SystemExit('unsupported file type: '+rp)
    with ProcessPoolExecutor(n) as ex:
        for e,h in zip(files, ex.map(sha,[os.path.join(root,e['path']) for e in files],chunksize=8)): e['sha256']=h
    files.sort(key=lambda e:e['path']); dirs.sort(key=lambda e:e['path'])
    json.dump({'root':root,'files':files,'dirs':dirs,'links':links,'bytes':sum(e['size'] for e in files)},open(out,'w'),indent=0)
    print(f'{root}: {len(files)} files, {sum(e["size"] for e in files)/1e9:.3f} GB')
if __name__=='__main__': main(sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 16)
