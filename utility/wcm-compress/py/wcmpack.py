#!/usr/bin/env python3
"""wcmpack — lossless archive codec for 4DWCM run directories (utility/wcm-compress).

    wcmpack.py pack   <run_dir>  <archive_dir> [--jobs N]
    wcmpack.py unpack <archive_dir> <out_dir>  [--jobs N]

Contract: unpack(pack(run)) reproduces every file byte for byte (plus mode and mtime), every directory and symlink.
Handlers claim files in registry order; whatever no handler claims goes through `raw` (concatenated, zstd).
Each handler must verify its own reconstruction at pack time and fall back to literal bytes where it cannot.

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)
"""
import os, sys, json, stat, time, hashlib, subprocess, argparse, mmap
from concurrent.futures import ProcessPoolExecutor

FORMAT = 1
JOBS = 24

# ---------------------------------------------------------------- streams
class StreamWriter:
    """Named append-only byte streams; each is compressed with its own command at close."""
    def __init__(self, arch):
        self.arch = arch; self.f = {}; self.cmd = {}; self.size = {}
    def open(self, name, cmd=('zstd', '-19', '--long=31')):
        if name not in self.f:
            self.f[name] = open(os.path.join(self.arch, name + '.raw'), 'wb'); self.cmd[name] = list(cmd); self.size[name] = 0
        return name
    def put(self, name, b):
        off = self.size[name]; self.f[name].write(b); self.size[name] += len(b); return off
    def close(self, jobs):
        info = {}
        for name, fh in self.f.items():
            fh.close(); raw = os.path.join(self.arch, name + '.raw'); cmd = self.cmd[name]
            t = time.time()
            if cmd[0] == 'zstd':
                subprocess.run(cmd + ['-T%d' % jobs, '-q', '-f', raw, '-o', raw[:-4] + '.zst'], check=True); ext = '.zst'
            elif cmd[0] == 'xz':
                subprocess.run(cmd + ['-T%d' % jobs, '-k', '-f', '-c', raw], check=True, stdout=open(raw[:-4] + '.xz', 'wb')); ext = '.xz'
            else:
                raise ValueError(cmd)
            os.remove(raw)
            info[name] = {'cmd': cmd, 'file': name + ext, 'raw_bytes': self.size[name],
                          'packed_bytes': os.path.getsize(os.path.join(self.arch, name + ext)), 'seconds': round(time.time() - t, 1)}
        return info

class StreamReader:
    def __init__(self, arch, info, scratch):
        self.arch = arch; self.info = info; self.scratch = scratch; self.mm = {}
    def get(self, name):
        if name not in self.mm:
            i = self.info[name]; src = os.path.join(self.arch, i['file']); dst = os.path.join(self.scratch, name + '.raw')
            if i['file'].endswith('.zst'):
                subprocess.run(['zstd', '-d', '-q', '-f', '--long=31', src, '-o', dst], check=True)
            else:
                subprocess.run(['xz', '-d', '-c', '-T%d' % JOBS, src], check=True, stdout=open(dst, 'wb'))
            fh = open(dst, 'rb')
            self.mm[name] = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) if os.path.getsize(dst) else b''
        return self.mm[name]
    def path(self, name):
        self.get(name); return os.path.join(self.scratch, name + '.raw')

# ---------------------------------------------------------------- handlers
class Raw:
    """Fallback: every unclaimed file's bytes, concatenated in path order, one zstd stream."""
    name = 'raw'
    def claim(self, root, entries): return list(entries)
    def pack(self, root, entries, sw, jobs):
        s = sw.open('misc')
        for e in entries:
            with open(os.path.join(root, e['path']), 'rb') as f:
                e['off'] = sw.put(s, f.read())
        return {}
    def unpack(self, out, entries, meta, sr, jobs):
        mm = sr.get('misc')
        for e in entries:
            write_file(out, e, mm[e['off']:e['off'] + e['size']])

HANDLERS = []          # specialised handlers, registry order; Raw is appended last

# ---------------------------------------------------------------- counts_fluxes_temp regenerated from counts_and_fluxes.csv
class CountsTemp:
    """counts_fluxes_temp/counts_fluxes_<t>.csv is column <t> of counts_and_fluxes.csv ("name,value" lines).
    claim() regenerates each file and keeps only those that match byte for byte; the rest fall through to raw."""
    name = 'cftemp'
    def _table(self, csv_path):
        rows = open(csv_path, 'rb').read().split(b'\n')
        if rows and rows[-1] == b'': rows = rows[:-1]
        return [r.split(b',') for r in rows]
    @staticmethod
    def render(tab, c):
        return b''.join(r[0] + b',' + r[c] + b'\n' for r in tab)
    def claim(self, root, entries):
        final = os.path.join(root, 'counts_and_fluxes.csv')
        cands = [e for e in entries if e['path'].startswith('counts_fluxes_temp/counts_fluxes_') and e['path'].endswith('.csv')]
        if not cands or not os.path.exists(final): return []
        tab = self._table(final); ncol = len(tab[0])
        if any(len(r) != ncol for r in tab): return []
        col = {h: i for i, h in enumerate(tab[0])}; mine = []
        for e in cands:
            t = e['path'][len('counts_fluxes_temp/counts_fluxes_'):-4].encode()
            if t in col and self.render(tab, col[t]) == open(os.path.join(root, e['path']), 'rb').read():
                e['col'] = t.decode(); mine.append(e)
        return mine
    def pack(self, root, entries, sw, jobs): return {'source': 'counts_and_fluxes.csv'}
    def unpack(self, out, entries, meta, sr, jobs):
        tab = self._table(os.path.join(out, meta['source'])); col = {h: i for i, h in enumerate(tab[0])}
        for e in entries: write_file(out, e, self.render(tab, col[e['col'].encode()]))

HANDLERS.append(CountsTemp())

# ---------------------------------------------------------------- float64 monomer coordinates, planar + byte-shuffled
_MONO_RE = __import__('re').compile(r'^DNA/dna_monomers_(\d+)\.bin$')

def _mono_steps(paths):
    return sorted(int(m.group(1)) for m in (_MONO_RE.match(p) for p in paths) if m)

def _shuffle(b):
    """per plane, the IEEE bit patterns mapped to order-preserving uint64 keys, second differences along the chain
    (mod 2^64, exactly invertible), zigzag, bytes of equal significance together."""
    import numpy as np
    u = np.ascontiguousarray(np.frombuffer(b, dtype='<u8').reshape(-1, 3).T)
    k = np.where(u >> np.uint64(63), ~u, u | np.uint64(1 << 63))
    with np.errstate(over='ignore'):
        r = np.diff(np.diff(k, axis=1, prepend=np.uint64(0)), axis=1, prepend=np.uint64(0)).view(np.int64)
    z = ((r << 1) ^ (r >> 63)).view(np.uint64)
    return np.ascontiguousarray(z).view(np.uint8).reshape(3, -1, 8).transpose(0, 2, 1).tobytes()

def _unshuffle(b):
    import numpy as np
    n = len(b) // 24
    z = np.frombuffer(b, dtype=np.uint8).reshape(3, 8, n).transpose(0, 2, 1).copy().view('<u8').reshape(3, n)
    r = (z >> np.uint64(1)) ^ (np.uint64(0) - (z & np.uint64(1)))
    with np.errstate(over='ignore'):
        k = np.cumsum(np.cumsum(r, axis=1, dtype=np.uint64), axis=1, dtype=np.uint64)
    u = np.where(k >> np.uint64(63), k & np.uint64((1 << 63) - 1), ~k)
    return np.ascontiguousarray(u.T).tobytes()

class Mono:
    """DNA/dna_monomers_<step>.bin (n×3 float64): x/y/z planes, chain second differences of the bit patterns, one zstd stream.
    Full precision — every mantissa bit is kept."""
    name = 'mono'
    def claim(self, root, entries):
        return [e for e in entries if _MONO_RE.match(e['path']) and e['size'] % 24 == 0 and e['size'] > 0]
    def pack(self, root, entries, sw, jobs):
        s = sw.open('mono')
        for e in entries:
            b = open(os.path.join(root, e['path']), 'rb').read(); t = _shuffle(b)
            assert _unshuffle(t) == b
            e['off'] = sw.put(s, t)
        return {}
    def unpack(self, out, entries, meta, sr, jobs):
        mm = sr.get('mono')
        for e in entries: write_file(out, e, _unshuffle(mm[e['off']:e['off'] + e['size']]))

HANDLERS.append(Mono())

# ---------------------------------------------------------------- HDF5 deflate chunks re-derived (+ sparse particle lattice)
H5SIG = b'\x89HDF\r\n\x1a\n'

def _h5_base(path):
    with open(path, 'rb') as f:
        for b in [0] + [512 << k for k in range(16)]:
            f.seek(b); sig = f.read(8)
            if len(sig) < 8: return None
            if sig == H5SIG: return b
    return None

def _h5_chunks(m, base):
    """All raw-data chunks referenced by level-0 v1 B-tree ('TREE', type 1) nodes: (addr, size). Key width is found per node by
    requiring the first child to be one complete zlib stream; overlapping hits are dropped (they stay in the skeleton)."""
    import struct, zlib
    N = len(m); pos = 0; ch = {}
    while True:
        i = m.find(b'TREE', pos)
        if i < 0: break
        pos = i + 4
        if i + 24 > N or m[i + 4] != 1 or m[i + 5] != 0: continue
        nent = struct.unpack_from('<H', m, i + 6)[0]
        if not 0 < nent <= 4096: continue
        for rank1 in (2, 3, 4, 5, 6, 7, 8):
            ks = 8 + 8 * rank1; ent = []; ok = True
            for e in range(nent):
                k = i + 24 + e * (ks + 8)
                if k + ks + 8 > N: ok = False; break
                csz = struct.unpack_from('<I', m, k)[0]; addr = struct.unpack_from('<Q', m, k + ks)[0] + base
                if struct.unpack_from('<Q', m, k + ks - 8)[0] != 0 or not (0 < csz and addr + csz <= N): ok = False; break
                ent.append((addr, csz))
            if ok:
                try:
                    d = zlib.decompressobj(); d.decompress(m[ent[0][0]:ent[0][0] + ent[0][1]]); ok = d.eof and not d.unused_data
                except Exception: ok = False
            if ok:
                for a, c in ent: ch[a] = c
                break
    out = []; end = 0
    for a in sorted(ch):
        c = ch[a]
        if a < end:
            if out and out[-1][0] + out[-1][1] > a: out.pop()
            continue
        out.append((a, c)); end = a + c
    return out

LAT_BYTES = 262144                                     # uint32 [16,16,16,16] chunk of the particle lattice (4096 sites × 16 slots)

def _lat_perm(order):
    """species → rank by frequency (the stored order first, every other uint16 value after it, ascending)."""
    import numpy as np
    seen = np.zeros(65536, bool); seen[order] = True
    perm = np.concatenate([np.asarray(order, np.int64), np.flatnonzero(~seen)]).astype('<u2')
    rank = np.empty(65536, '<u2'); rank[perm] = np.arange(65536, dtype='<u2'); return rank, perm

def _lat_enc(d, rank):
    import numpy as np
    a = np.frombuffer(d, dtype='<u4').reshape(4096, 16); nz = a != 0; occ = nz.sum(1).astype(np.uint8)
    if not (nz == (np.arange(16)[None, :] < occ[:, None])).all() or a.max() > 0xFFFF: return None
    sp = rank[a[nz]]; return occ.tobytes(), sp.view(np.uint8)[1::2].tobytes(), sp.view(np.uint8)[0::2].tobytes()

def _lat_dec(occ, hi, lo, perm):
    import numpy as np
    o = np.frombuffer(occ, dtype=np.uint8); n = len(hi)
    sp = np.empty(2 * n, np.uint8); sp[0::2] = np.frombuffer(lo, np.uint8); sp[1::2] = np.frombuffer(hi, np.uint8)
    a = np.zeros((4096, 16), '<u4'); a[np.arange(16)[None, :] < o[:, None]] = perm[sp.view('<u2')]; return a.tobytes()

def _h5_count_job(args):
    import zlib, numpy as np
    path, chunks = args; fh = open(path, 'rb'); m = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ); cnt = np.zeros(65536, np.int64)
    for a, c in chunks:
        try: d = zlib.decompress(m[a:a + c])
        except Exception: continue
        if len(d) == LAT_BYTES:
            v = np.frombuffer(d, '<u4'); v = v[(v != 0) & (v <= 0xFFFF)]; cnt += np.bincount(v, minlength=65536)
    m.close(); return cnt

def _h5_pack_job(args):
    import zlib
    path, chunks, order = args; fh = open(path, 'rb'); m = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
    rank, perm = _lat_perm(order); res = []                                           # per chunk: (level or 0, kind, occ, hi, lo, raw)
    for a, c in chunks:
        raw = m[a:a + c]; r = (0, 0, b'', b'', b'', b'')
        try:
            o = zlib.decompressobj(); d = o.decompress(raw)
            if o.eof and not o.unused_data:
                lvl = next((l for l in (1, 6, 9, 2, 3, 4, 5, 7, 8, 0) if zlib.compress(d, l) == raw), None)
                if lvl is not None:
                    e = _lat_enc(d, rank) if len(d) == LAT_BYTES else None
                    if e is not None and _lat_dec(*e, perm) == d: r = (lvl, 1, e[0], e[1], e[2], b'')
                    else: r = (lvl, 2, b'', b'', b'', d)
        except Exception: pass
        res.append(r)
    m.close(); return res

def _h5_unpack_job(args):
    import zlib
    paths, recs, order = args; mm = {k: (open(p, 'rb') if p else None) for k, p in paths.items()}; perm = _lat_perm(order)[1]
    def rd(k, off, n):
        f = mm[k]; f.seek(off); return f.read(n)
    out = []
    for addr, lvl, kind, oo, no, so, ns, ro, nr in recs:
        d = _lat_dec(rd('occ', oo, no), rd('hi', so, ns), rd('lo', so, ns), perm) if kind == 1 else rd('raw', ro, nr)
        out.append((addr, zlib.compress(d, lvl)))
    return out

class H5:
    """HDF5 files (LM .lm, with or without a user block): every deflate chunk whose zlib re-compression reproduces its stored
    bytes is replaced by its content — particle-lattice chunks as (per-site occupancy, species-rank hi bytes, species-rank lo bytes),
    others verbatim — and re-deflated at unpack. The rest of the file (superblock, B-trees, headers, uncompressed datasets,
    non-reproducible chunks) is kept byte for byte in a skeleton stream."""
    name = 'h5'
    def claim(self, root, entries):
        return [e for e in entries if e['size'] > 1 << 20 and _h5_base(os.path.join(root, e['path'])) is not None]
    def pack(self, root, entries, sw, jobs):
        import numpy as np
        S = {k: sw.open('h5_' + k) for k in ('skel', 'occ', 'hi', 'lo', 'raw', 'index')}; stats = {}
        with ProcessPoolExecutor(jobs) as ex:
            for e in entries:
                p = os.path.join(root, e['path']); base = _h5_base(p); fh = open(p, 'rb'); m = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
                ch = _h5_chunks(m, base); step = max(1, len(ch) // (jobs * 16) + 1)
                parts = [ch[i:i + step] for i in range(0, len(ch), step)]
                import numpy as np
                cnt = sum(ex.map(_h5_count_job, [(p, q) for q in parts[::8]]), np.zeros(65536, np.int64)) if parts else np.zeros(65536, np.int64)
                order = [int(v) for v in np.argsort(-cnt, kind='stable') if cnt[v] > 0]; e['rank'] = order
                keep = []; lv = []; kd = []; no = []; ns = []; nr = []
                for part, res in zip(parts, ex.map(_h5_pack_job, [(p, q, order) for q in parts])):
                    for (a, c), (l, k, occ, hi, lo, raw) in zip(part, res):
                        if not k: continue
                        keep.append((a, c)); lv.append(l); kd.append(k); no.append(len(occ)); ns.append(len(hi)); nr.append(len(raw))
                        if k == 1: sw.put(S['occ'], occ); sw.put(S['hi'], hi); sw.put(S['lo'], lo)
                        else: sw.put(S['raw'], raw)
                cur = 0; e['skel'] = sw.size[S['skel']]
                for a, c in keep: sw.put(S['skel'], m[cur:a]); cur = a + c
                sw.put(S['skel'], m[cur:]); m.close()
                e['base'] = [sw.size[S[k]] for k in ('occ', 'hi', 'raw')]
                e['base'][0] -= sum(no); e['base'][1] -= sum(ns); e['base'][2] -= sum(nr)
                idx = np.array(keep, dtype='<u8').reshape(-1, 2)
                cols = [idx[:, 0], idx[:, 1], np.array(lv, '<u8'), np.array(kd, '<u8'), np.array(no, '<u8'), np.array(ns, '<u8'), np.array(nr, '<u8')]
                blob = np.stack(cols).astype('<u8').tobytes() if keep else b''
                e['index'] = (sw.put(S['index'], blob), len(keep))
                stats[e['path']] = {'chunks_found': len(ch), 'chunks_rederived': len(keep), 'lattice': int(sum(1 for k in kd if k == 1))}
        return stats
    def unpack(self, out, entries, meta, sr, jobs):
        import numpy as np
        paths = {k: sr.path('h5_' + v) if v in sr.info or ('h5_' + v) in sr.info else None for k, v in (('occ', 'occ'), ('hi', 'hi'), ('lo', 'lo'), ('raw', 'raw'))}
        skel = sr.get('h5_skel'); ix = sr.get('h5_index')
        with ProcessPoolExecutor(jobs) as ex:
            for e in entries:
                off, n = e['index']; A = np.frombuffer(ix[off:off + 56 * n], '<u8').reshape(7, n) if n else np.zeros((7, 0), '<u8')
                addr, csz, lvl, kind, no, ns, nr = A
                oo = e['base'][0] + np.concatenate([[0], np.cumsum(no)[:-1]]) if n else []
                so = e['base'][1] + np.concatenate([[0], np.cumsum(ns)[:-1]]) if n else []
                ro = e['base'][2] + np.concatenate([[0], np.cumsum(nr)[:-1]]) if n else []
                p = os.path.join(out, e['path']); os.makedirs(os.path.dirname(p), exist_ok=True)
                with open(p, 'wb') as f:
                    f.truncate(e['size']); cur = 0; sk = e['skel']
                    for a, c in zip(addr.tolist(), csz.tolist()):
                        f.seek(cur); f.write(skel[sk:sk + a - cur]); sk += a - cur; cur = a + c
                    f.seek(cur); f.write(skel[sk:sk + e['size'] - cur])
                    recs = list(zip(addr.tolist(), lvl.tolist(), kind.tolist(), list(map(int, oo)), no.tolist(), list(map(int, so)), ns.tolist(), list(map(int, ro)), nr.tolist()))
                    step = max(1, len(recs) // (jobs * 16) + 1)
                    for res in ex.map(_h5_unpack_job, [(paths, recs[i:i + step], e.get('rank', [])) for i in range(0, len(recs), step)]):
                        for a, b in res:
                            f.seek(a); f.write(b)

HANDLERS.append(H5())

# ---------------------------------------------------------------- DNA text files regenerated from dna_monomers_<step>.bin
def _load_mono(root, step):
    import numpy as np
    return np.fromfile(os.path.join(root, 'DNA/dna_monomers_%d.bin' % step), dtype='<f8').reshape(-1, 3).tolist()

def _enc_data_lines(lines, mono):
    """data.lammps Atoms lines 'id  \tmol\ttype\tx  \ty  \tz' whose x y z == '%g' of mono[id-1] → 'id  \tmol\ttype\x00'."""
    out = []; hit = 0; inatoms = False; seen = False
    for l in lines:
        if l.startswith('Atoms'): inatoms = True; seen = False; out.append(l); continue
        if inatoms:
            if l == '':
                if seen: inatoms = False
                out.append(l); continue
            seen = True
            f = l.split('\t')
            if len(f) == 6 and mono is not None:
                try: i = int(f[0]) - 1
                except ValueError: i = -1
                if 0 <= i < len(mono):
                    x = mono[i]
                    if l == '%s\t%s\t%s\t%g  \t%g  \t%g' % (f[0], f[1], f[2], x[0], x[1], x[2]):
                        out.append('%s\t%s\t%s\x00' % (f[0], f[1], f[2])); hit += 1; continue
        out.append(l)
    return out, hit

def _dec_data_lines(lines, mono):
    out = []
    for l in lines:
        if l.endswith('\x00'):
            p = l[:-1]; x = mono[int(p.split('\t', 1)[0]) - 1]; out.append('%s\t%g  \t%g  \t%g' % (p, x[0], x[1], x[2]))
        else: out.append(l)
    return out

def _enc_traj_lines(lines, mono):
    """lammpstrj atom lines 'id type x y z cid ctype' with x y z == '%g' of mono[id-1] → 'id type\x00 cid ctype'."""
    out = []; hit = 0
    for l in lines:
        f = l.split(' ')
        if len(f) == 7 and mono is not None and not l.startswith('ITEM'):
            try: i = int(f[0]) - 1
            except ValueError: i = -1
            if 0 <= i < len(mono):
                x = mono[i]
                if l == '%s %s %g %g %g %s %s' % (f[0], f[1], x[0], x[1], x[2], f[5], f[6]):
                    out.append('%s %s\x00%s %s' % (f[0], f[1], f[5], f[6])); hit += 1; continue
        out.append(l)
    return out, hit

def _dec_traj_lines(lines, mono):
    out = []
    for l in lines:
        if '\x00' in l:
            a, b = l.split('\x00'); x = mono[int(a.split(' ', 1)[0]) - 1]; out.append('%s %g %g %g %s' % (a, x[0], x[1], x[2], b))
        else: out.append(l)
    return out

def _dna_data_job(args):
    root, rel, cands = args
    txt = open(os.path.join(root, rel), 'rb').read()
    if b'\x00' in txt: return None
    lines = txt.decode('latin-1').split('\n'); best = (lines, 0, None)
    for st in cands:                                   # the step before first (that is what btree wrote), then the same step
        enc, hit = _enc_data_lines(lines, _load_mono(root, st))
        if hit > best[1]: best = (enc, hit, st)
        if hit: break
    enc, hit, st = best
    if st is not None and '\n'.join(_dec_data_lines(enc, _load_mono(root, st))).encode('latin-1') != txt: return None
    return '\n'.join(enc).encode('latin-1'), st, hit

def _dna_traj_job(args):
    root, rel, a, b, cands = args
    with open(os.path.join(root, rel), 'rb') as f:
        f.seek(a); txt = f.read(b - a)
    lines = txt.decode('latin-1').split('\n'); best = (lines, 0, None)
    if b'\x00' not in txt:
        for st in cands:
            if st is None: continue
            enc, hit = _enc_traj_lines(lines, _load_mono(root, st))
            if hit > best[1]: best = (enc, hit, st)
            if hit: break
    enc, hit, st = best
    if st is not None and '\n'.join(_dec_traj_lines(enc, _load_mono(root, st))).encode('latin-1') != txt: enc, st = lines, None
    return '\n'.join(enc).encode('latin-1'), st

def _dna_data_unjob(args):
    out, skel_path, off, n, st = args
    with open(skel_path, 'rb') as f:
        f.seek(off); sk = f.read(n)
    if st is None: return sk
    return '\n'.join(_dec_data_lines(sk.decode('latin-1').split('\n'), _load_mono(out, st))).encode('latin-1')

def _dna_traj_unjob(args):
    out, skel_path, off, n, st = args
    with open(skel_path, 'rb') as f:
        f.seek(off); sk = f.read(n)
    if st is None: return sk
    return '\n'.join(_dec_traj_lines(sk.decode('latin-1').split('\n'), _load_mono(out, st))).encode('latin-1')

class DnaText:
    """DNA/data.lammps_<N>: DNA coordinates == '%g' of dna_monomers_<prev step>.bin; DNA/chromosome.lammpstrj frame k: DNA
    coordinates == '%g' of dna_monomers_<steps[k-1]>.bin. The coordinates are dropped and rebuilt; the remaining skeletons
    (ids, types, boundary shell, bonds, angles, frame headers) are de-duplicated by content and zstd'd in one long-window stream.
    Every file / frame is decoded again at pack time and compared with the original; a mismatch keeps it literal."""
    name = 'dnatext'
    def claim(self, root, entries):
        d = os.path.join(root, 'DNA'); self.steps = _mono_steps('DNA/' + f for f in (os.listdir(d) if os.path.isdir(d) else []))
        if not self.steps: return []
        return [e for e in entries if __import__('re').match(r'^DNA/data\.lammps_\d+$', e['path']) or e['path'] == 'DNA/chromosome.lammpstrj']
    def pack(self, root, entries, sw, jobs):
        import bisect
        steps = self.steps; s = sw.open('dna_skel'); seen = {}
        def put(b):
            h = hashlib.sha256(b).hexdigest()
            if h not in seen: seen[h] = (sw.put(s, b), len(b))
            return seen[h]
        stats = {'data_files': 0, 'data_regen': 0, 'frames': 0, 'frames_regen': 0, 'lines': 0}
        data = [e for e in entries if e['path'] != 'DNA/chromosome.lammpstrj']; traj = [e for e in entries if e['path'] == 'DNA/chromosome.lammpstrj']
        jobsl = []
        for e in data:
            n = int(e['path'].rsplit('_', 1)[1]); k = bisect.bisect_left(steps, n)
            jobsl.append((root, e['path'], [x for x in (steps[k - 1] if k > 0 else None, n if k < len(steps) and steps[k] == n else None) if x is not None]))
        with ProcessPoolExecutor(jobs) as ex:
            for e, r in zip(data, ex.map(_dna_data_job, jobsl, chunksize=4)):
                if r is None:
                    e['skel'] = put(open(os.path.join(root, e['path']), 'rb').read()); e['src'] = None
                else:
                    e['skel'] = put(r[0]); e['src'] = r[1]; stats['data_regen'] += r[1] is not None; stats['lines'] += r[2]
                stats['data_files'] += 1
            for e in traj:
                p = os.path.join(root, e['path']); fh = open(p, 'rb'); m = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
                offs = []; q = 0
                while True:
                    i = m.find(b'ITEM: TIMESTEP\n', q)
                    if i < 0: break
                    offs.append(i); q = i + 1
                if not offs or offs[0] != 0: offs = [0] + offs
                bounds = list(zip(offs, offs[1:] + [len(m)])); m.close()
                tj = [(root, e['path'], a, b, [steps[k - 1] if 0 < k <= len(steps) else None]) for k, (a, b) in enumerate(bounds)]
                e['frames'] = []
                for r in ex.map(_dna_traj_job, tj, chunksize=2):
                    e['frames'].append((put(r[0]), r[1])); stats['frames'] += 1; stats['frames_regen'] += r[1] is not None
        stats['unique_skeletons'] = len(seen)
        return stats
    def unpack(self, out, entries, meta, sr, jobs):
        sp = sr.path('dna_skel')
        with ProcessPoolExecutor(jobs) as ex:
            data = [e for e in entries if 'frames' not in e]
            for e, b in zip(data, ex.map(_dna_data_unjob, [(out, sp, e['skel'][0], e['skel'][1], e['src']) for e in data], chunksize=4)):
                write_file(out, e, b)
            for e in entries:
                if 'frames' not in e: continue
                p = os.path.join(out, e['path']); os.makedirs(os.path.dirname(p), exist_ok=True)
                with open(p, 'wb') as f:
                    for b in ex.map(_dna_traj_unjob, [(out, sp, o, n, st) for (o, n), st in e['frames']], chunksize=2): f.write(b)

HANDLERS.append(DnaText())

# ---------------------------------------------------------------- tree walk / metadata
def write_file(out, e, data):
    p = os.path.join(out, e['path']); os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, 'wb') as f: f.write(data)

def walk(root):
    files, dirs, links = [], [], []
    for d, ds, fs in os.walk(root):
        rd = os.path.relpath(d, root); st = os.lstat(d)
        dirs.append({'path': rd, 'mode': stat.S_IMODE(st.st_mode), 'mtime_ns': st.st_mtime_ns})
        for f in fs:
            p = os.path.join(d, f); st = os.lstat(p); rp = os.path.relpath(p, root)
            if stat.S_ISLNK(st.st_mode): links.append({'path': rp, 'target': os.readlink(p)})
            elif stat.S_ISREG(st.st_mode): files.append({'path': rp, 'size': st.st_size, 'mode': stat.S_IMODE(st.st_mode), 'mtime_ns': st.st_mtime_ns})
            else: raise SystemExit('unsupported file type: ' + rp)
    files.sort(key=lambda e: e['path']); dirs.sort(key=lambda e: e['path'])
    return files, dirs, links

def _sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''): h.update(b)
    return h.hexdigest()

# ---------------------------------------------------------------- pack / unpack
def pack(root, arch, jobs):
    t0 = time.time(); root = os.path.abspath(root); os.makedirs(arch, exist_ok=True)
    files, dirs, links = walk(root)
    with ProcessPoolExecutor(jobs) as ex:
        for e, h in zip(files, ex.map(_sha, [os.path.join(root, e['path']) for e in files], chunksize=8)): e['sha256'] = h
    sw = StreamWriter(arch); left = {e['path']: e for e in files}; hmeta = {}; timing = {}
    for h in HANDLERS + [Raw()]:
        t = time.time()
        mine = h.claim(root, [left[p] for p in sorted(left)])
        for e in mine: e['h'] = h.name; del left[e['path']]
        hmeta[h.name] = h.pack(root, mine, sw, jobs) if mine else {}
        timing[h.name] = round(time.time() - t, 1)
    t = time.time(); sinfo = sw.close(jobs); timing['compress'] = round(time.time() - t, 1)
    man = {'format': FORMAT, 'source': root, 'files': files, 'dirs': dirs, 'links': links, 'handlers': hmeta,
           'streams': sinfo, 'order': [h.name for h in [Raw()] + HANDLERS], 'timing': timing}
    raw = json.dumps(man).encode()
    subprocess.run(['zstd', '-19', '-q', '-f', '-o', os.path.join(arch, 'manifest.json.zst')], input=raw, check=True)
    total = sum(os.path.getsize(os.path.join(arch, f)) for f in os.listdir(arch))
    print(json.dumps({'pack_s': round(time.time() - t0, 1), 'raw_bytes': sum(e['size'] for e in files), 'packed_bytes': total,
                      'timing': timing, 'streams': {k: (v['raw_bytes'], v['packed_bytes']) for k, v in sinfo.items()}}))

def unpack(arch, out, jobs):
    t0 = time.time()
    man = json.loads(subprocess.run(['zstd', '-d', '-q', '-c', os.path.join(arch, 'manifest.json.zst')], capture_output=True, check=True).stdout)
    assert man['format'] == FORMAT
    os.makedirs(out, exist_ok=True); scratch = os.path.join(out, '.wcmpack_tmp'); os.makedirs(scratch, exist_ok=True)
    sr = StreamReader(arch, man['streams'], scratch)
    for d in man['dirs']: os.makedirs(os.path.join(out, d['path']), exist_ok=True)
    reg = {h.name: h for h in HANDLERS + [Raw()]}
    for name in man['order']:                       # raw first: derived handlers may read restored files
        es = [e for e in man['files'] if e.get('h') == name]
        if es: reg[name].unpack(out, es, man['handlers'][name], sr, jobs)
    for m in sr.mm.values():
        if m: m.close()
    for f in os.listdir(scratch): os.remove(os.path.join(scratch, f))
    os.rmdir(scratch)
    for l in man['links']: os.symlink(l['target'], os.path.join(out, l['path']))
    for e in man['files']:
        p = os.path.join(out, e['path']); os.chmod(p, e['mode']); os.utime(p, ns=(e['mtime_ns'], e['mtime_ns']))
    for d in sorted(man['dirs'], key=lambda d: -d['path'].count('/')):
        p = os.path.join(out, d['path']); os.chmod(p, d['mode']); os.utime(p, ns=(d['mtime_ns'], d['mtime_ns']))
    print(json.dumps({'unpack_s': round(time.time() - t0, 1)}))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['pack', 'unpack']); ap.add_argument('src'); ap.add_argument('dst')
    ap.add_argument('--jobs', type=int, default=JOBS); a = ap.parse_args(); JOBS = a.jobs
    (pack if a.cmd == 'pack' else unpack)(a.src, a.dst, a.jobs)
