# cython: boundscheck=False, wraparound=False, cdivision=True, language_level=3
"""the ribosome hook's common path (no relocation) as compiled passes over the lattice masks.

Authors
-------
Ron Acda — using an iterative LLM-guided workflow
    (https://github.com/quarkron/iterative-hillclimber/tree/main)

Same results as RibosomesRDME.placeRibosomes + updateRiboSites (checked with WCM_RIBO_COMPILED_VERIFY=1):
  coords2(a, b)            -> np.vstack([_true_coords(a), _true_coords(b)]) (rows in C order of each mask, a first)
  new_masks(pts, shape)    -> the new centre mask (True at pts) and cross mask (the six clipped neighbours of every centre)
  site_writes(...)         -> the final site types updateRiboSites leaves in the site-lattice view: its write groups in order
                              (old centres by cyto > outer > DNA, old crosses by cyto > outer > DNA > extracellular, new crosses off
                              membrane and DNA, new centres), later groups overwriting earlier ones, and the last site the Python
                              code wrote (for its single setSiteType call).
Masks are C-contiguous bool arrays read as flat bytes; 8 sites are tested at a time as one 64-bit word (the masks are sparse).
"""
import numpy as np
cimport numpy as cnp
from libc.stdint cimport uint8_t, uint64_t
from libc.string cimport memcpy


cdef inline uint64_t w8(const uint8_t* p, Py_ssize_t i) nogil:
    cdef uint64_t w
    memcpy(&w, p + i, 8)
    return w


cdef Py_ssize_t collect(const uint8_t* a, Py_ssize_t n, Py_ssize_t* out) nogil:
    """flat indices of the nonzero bytes of a[0:n], ascending"""
    cdef Py_ssize_t i = 0, r = 0, j
    while i + 8 <= n:
        if w8(a, i) != 0:
            for j in range(i, i + 8):
                if a[j]:
                    out[r] = j; r += 1
        i += 8
    while i < n:
        if a[i]:
            out[r] = i; r += 1
        i += 1
    return r


def coords2(cnp.ndarray a, cnp.ndarray b):
    cdef Py_ssize_t n = a.size, n1 = a.shape[1], n2 = a.shape[2], ra, rb, r
    ia = np.empty(n, dtype=np.intp); ib = np.empty(n, dtype=np.intp)
    cdef cnp.intp_t[::1] va = ia, vb = ib
    ra = collect(<const uint8_t*>cnp.PyArray_DATA(a), n, <Py_ssize_t*>&va[0])
    rb = collect(<const uint8_t*>cnp.PyArray_DATA(b), n, <Py_ssize_t*>&vb[0])
    out = np.empty((ra + rb, 3), dtype=np.intp)
    cdef cnp.intp_t[:, ::1] o = out
    cdef Py_ssize_t q, f, n12 = n1 * n2
    for r in range(ra):
        f = va[r]; q = f // n12; o[r, 0] = q; f -= q * n12; o[r, 1] = f // n2; o[r, 2] = f - (f // n2) * n2
    for r in range(rb):
        f = vb[r]; q = f // n12; o[ra + r, 0] = q; f -= q * n12; o[ra + r, 1] = f // n2; o[ra + r, 2] = f - (f // n2) * n2
    return out, ra


def new_masks(cnp.intp_t[:, ::1] pts, Py_ssize_t n0, Py_ssize_t n1, Py_ssize_t n2):
    centers = np.zeros((n0, n1, n2), dtype=np.bool_)
    crosses = np.zeros((n0, n1, n2), dtype=np.bool_)
    cdef uint8_t[:, :, ::1] c = centers.view(np.uint8), x = crosses.view(np.uint8)
    cdef Py_ssize_t p, i, j, k
    for p in range(pts.shape[0]):
        i = pts[p, 0]; j = pts[p, 1]; k = pts[p, 2]
        c[i, j, k] = 1
        x[i + 1 if i + 1 < n0 else n0 - 1, j, k] = 1
        x[i - 1 if i > 0 else 0, j, k] = 1
        x[i, j + 1 if j + 1 < n1 else n1 - 1, k] = 1
        x[i, j - 1 if j > 0 else 0, k] = 1
        x[i, j, k + 1 if k + 1 < n2 else n2 - 1] = 1
        x[i, j, k - 1 if k > 0 else 0] = 1
    return centers, crosses


def site_writes(cnp.ndarray svarr, cnp.ndarray old_c_, cnp.ndarray old_x_, cnp.ndarray new_c_, cnp.ndarray new_x_,
                cnp.ndarray cyto_, cnp.ndarray outer_, cnp.ndarray dna_, cnp.ndarray extra_, cnp.ndarray memb_,
                int i_cyto, int i_outer, int i_dna, int i_extra, int i_cross, int i_center):
    """Writes the final values into the site view; returns the (x, y, z) of the site the Python code wrote last, or None."""
    cdef const uint8_t* oc = <const uint8_t*>cnp.PyArray_DATA(old_c_)
    cdef const uint8_t* ox = <const uint8_t*>cnp.PyArray_DATA(old_x_)
    cdef const uint8_t* nc = <const uint8_t*>cnp.PyArray_DATA(new_c_)
    cdef const uint8_t* nx = <const uint8_t*>cnp.PyArray_DATA(new_x_)
    cdef const uint8_t* cy = <const uint8_t*>cnp.PyArray_DATA(cyto_)
    cdef const uint8_t* ou = <const uint8_t*>cnp.PyArray_DATA(outer_)
    cdef const uint8_t* dn = <const uint8_t*>cnp.PyArray_DATA(dna_)
    cdef const uint8_t* ex = <const uint8_t*>cnp.PyArray_DATA(extra_)
    cdef const uint8_t* me = <const uint8_t*>cnp.PyArray_DATA(memb_)
    cdef uint8_t[:, :, :] sv = svarr
    cdef Py_ssize_t n = old_c_.size, n1 = old_c_.shape[1], n2 = old_c_.shape[2], n12 = n1 * n2, i = 0, e, j, f
    # g0 old centres cyto, g1 outer, g2 DNA, g3 old crosses cyto, g4 outer, g5 DNA, g6 extracellular, g7 new crosses, g8 new centres
    cdef Py_ssize_t last[9]
    cdef int g, h, v
    for g in range(9): last[g] = -1
    while i < n:
        e = i + 8 if i + 8 <= n else n
        if e - i == 8 and (w8(oc, i) | w8(ox, i) | w8(nc, i) | w8(nx, i)) == 0:
            i = e; continue
        for j in range(i, e):
            v = -1
            if oc[j]:
                g = -1
                if cy[j]: v = i_cyto; g = 0
                elif ou[j]: v = i_outer; g = 1
                elif dn[j]: v = i_dna; g = 2
                if g >= 0: last[g] = j
            if ox[j]:
                g = -1
                if cy[j]: v = i_cyto; g = 3
                elif ou[j]: v = i_outer; g = 4
                elif dn[j]: v = i_dna; g = 5
                elif ex[j]: v = i_extra; g = 6
                if g >= 0: last[g] = j
            if nx[j] and not me[j] and not dn[j]:
                v = i_cross; last[7] = j
            if nc[j]:
                v = i_center; last[8] = j
            if v >= 0:
                f = j
                sv[f // n12, (f % n12) // n2, f % n2] = <uint8_t>v
        i = e
    for h in range(8, -1, -1):
        if last[h] >= 0:
            f = last[h]
            return (f // n12, (f % n12) // n2, f % n2)
    return None
