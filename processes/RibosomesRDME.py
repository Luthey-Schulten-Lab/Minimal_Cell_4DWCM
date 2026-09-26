"""
Update ribosome excluded-volume sites on the RDME lattice.

Authors
-------
Ron Acda — vectorised ribosome scans and placement, compiled common path (ribo_fast.pyx),
    GPU-side particle-lattice reads
    (using an iterative LLM-guided workflow: https://github.com/quarkron/iterative-hillclimber/tree/main)
Alfia Parvez — KDTree cache, vectorized center/cross updates, faster site walks
Zane Thornburg — original ribosome placement and site update routines
"""

import numpy as np
from scipy import spatial
from utility.LatticeFunctions import *

# KDTree over outer-cytoplasm shell; rebuild when shell point count changes.
_cached_tree = None
_cached_she_coords_size = None

def _coords(mask):
    """``np.argwhere(mask)`` for a 3-D mask: same rows, same (C) order, same intp dtype, ~3x faster on the
    64x64x128 region masks (flatnonzero + unravel_index instead of the generic N-d nonzero)."""
    return np.column_stack(np.unravel_index(np.flatnonzero(mask), mask.shape))


def _true_coords(mask):
    """``_coords(mask == True)`` without building the boolean copy when the mask already is boolean (same rows, same order)."""
    return _coords(mask if mask.dtype == np.bool_ else mask == True)


# compiled passes for the common (no-relocation) path; WCM_RIBO_COMPILED_OFF=1 = the numpy code, WCM_RIBO_COMPILED_VERIFY=1 runs both
import os as _os_env
_RIBO_COMPILED_OFF = _os_env.environ.get('WCM_RIBO_COMPILED_OFF') is not None
_RIBO_COMPILED_VERIFY = _os_env.environ.get('WCM_RIBO_COMPILED_VERIFY') is not None
_w127 = None
if not _RIBO_COMPILED_OFF:
    try:
        import pyximport as _pyx
        _pyx.install(setup_args={'include_dirs': np.get_include()}, language_level=3)
        from processes import ribo_fast as _w127
    except Exception as _e:   # no compiler / build failure: the numpy path, as before
        print('compiled ribosome path unavailable (%s), using numpy' % _e)
        _w127 = None


def _w127_ok(*masks):
    return _w127 is not None and all(m.dtype == np.bool_ and m.ndim == 3 and m.flags.c_contiguous for m in masks)


def _ribo_sites(plattice, ridx):
    """``(plattice == ridx).any(axis=0)[..., 0]``: sites holding particle type ``ridx`` in any slot.
    Lattice Microbes keeps each site's particles packed from slot 0 (getOccupancy counts up to the first empty slot; the RDME
    kernels, addParticle and deleteParticle all keep it), so slot s can only be occupied where slot s-1 is. Slot 0 is scanned
    in full and each further slot only at the sites still occupied: ~5 MB of the particle lattice read instead of all 33.5 MB
    (3.6 -> ~1 ms per hook from a cold cache), same mask."""
    s0 = plattice[0, ..., 0]
    m = s0 == ridx
    mf = m.reshape(-1)
    idx = np.flatnonzero(s0 != 0)   # the same indices; nonzero of a bool mask is ~3x faster than of uint32
    for s in range(1, plattice.shape[0]):
        if idx.size == 0:
            break
        v = plattice[s, ..., 0].reshape(-1)[idx]
        mf[idx[v == ridx]] = True
        idx = idx[v != 0]
    return m


#########################################################################################
wcm_last_placement_touched_particles = True   # set by every placeRibosomes call


def placeRibosomes(lattice, sim_properties, region_dict, ribo_site_dict, updateTranslat=True, growth_step=False):
    """
    Place / relocate ribosome excluded-volume centers and crosses on the lattice.
    """
    global _cached_tree, _cached_she_coords_size, wcm_last_placement_touched_particles
    # whether this call may have changed the particle lattice (relocations, or a growth-step call); the hook then
    # tells Lattice Microbes to upload the particles as before, otherwise only the site types
    wcm_last_placement_touched_particles = bool(growth_step)
    
    N_edges = sim_properties['lattice_edges']
    ribo_IDs = sim_properties['riboIDs']
    
    membrane = region_dict['membrane']['shape']
    extracellular = region_dict['extracellular']['shape']
    cyto_shell = region_dict['outer_cytoplasm']['shape']
    
    # when the solver deferred the particle download for this hook (host particles stale), the two reads below run
    # on the device (the same mask and slots) and the particle view, which downloads the particles, is taken only if a
    # ribosome has to be relocated. WCM_RIBO_DEVICE_VERIFY=1 also computes both from the host copy and compares.
    on_device = _wcm_device_reads(lattice)
    plattice = None if on_device else lattice.getParticleLatticeView()
    
    ribo_center_points = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    
    # Rebuild KDTree when outer-cytoplasm shell size changes.
    # the cache test only needs the shell's size; the coordinates themselves (a flatnonzero over the whole mask) are
    # built only where they are used: to rebuild the tree, and in the relocation fallback (from the current shell, as before).
    current_size = int(np.count_nonzero(cyto_shell if cyto_shell.dtype == np.bool_ else cyto_shell == True))
    she_coords = None
    
    if (_cached_tree is None or _cached_she_coords_size != current_size):
        if current_size > 0:
            she_coords = _true_coords(cyto_shell)
            tree = spatial.KDTree(she_coords)
            _cached_tree = tree
            _cached_she_coords_size = current_size
        else:
            # Empty shell — clear EV masks and return.
            ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
            ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
            return ribo_site_dict
    else:
        tree = _cached_tree
    
    RIBOidx = sim_properties['name_to_index']['ribosomeP']
    # Sites holding a ribosome particle in any slot: reduce over the slot axis first (one pass over the
    # particle lattice) instead of argwhere over every slot (16.8 -> ~1.5 ms per hook, same mask).
    if on_device:
        ribo_mask = _wcm_device_ribo_mask(lattice, RIBOidx, ribo_site_dict['ribos']['centers'].shape)
    else:
        ribo_mask = _ribo_sites(plattice, RIBOidx)   # packed-slot scan, same mask
    ribo_site_dict['ribos']['centers'] |= ribo_mask
    
    old_ribo_sites_list = []
    _t127 = next(iter(ribo_site_dict.values())) if len(ribo_site_dict) == 1 else None
    if _t127 is not None and _w127_ok(_t127['centers'], _t127['crosses']):
        _o, _nc = _w127.coords2(_t127['centers'], _t127['crosses'])
        if _RIBO_COMPILED_VERIFY:
            _ref = [a for a in (_true_coords(_t127['centers']), _true_coords(_t127['crosses'])) if len(a)]
            assert (np.vstack(_ref) if _ref else np.empty((0, 3), np.intp)).tolist() == _o.tolist(), 'compiled ribosome pass: coords2 mismatch'
        if len(_o):
            old_ribo_sites_list.append(_o)
        _t127 = 'done'
    for ribo_type, type_dict in (ribo_site_dict.items() if _t127 != 'done' else ()):
        ribo_center_sites = _true_coords(type_dict['centers'])
        if len(ribo_center_sites) > 0:
            old_ribo_sites_list.append(ribo_center_sites)
        
        ribo_cross_sites = _true_coords(type_dict['crosses'])
        if len(ribo_cross_sites) > 0:
            old_ribo_sites_list.append(ribo_cross_sites)
    
    if not old_ribo_sites_list:
        ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
        ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
        return ribo_site_dict
    
    old_ribo_sites = np.vstack(old_ribo_sites_list)
    
    # Fast path (every captured hook): no ribosome sits on a membrane / extracellular site, so nothing is
    # relocated and the result is just the old sites that hold a ribosome particle, found with one gather
    # of their slots. Any hook that needs a relocation runs the original per-site loop unchanged, because
    # relocations mutate the lattice while the loop is still reading it.
    xs, ys, zs = old_ribo_sites[:, 0], old_ribo_sites[:, 1], old_ribo_sites[:, 2]
    if on_device:
        slots = _wcm_device_slots(lattice, xs, ys, zs, ribo_site_dict['ribos']['centers'].shape)
        if _RIBO_DEVICE_VERIFY:
            _wcm_116_check(lattice, RIBOidx, ribo_mask, slots, xs, ys, zs)
    else:
        slots = plattice[:, xs, ys, zs, 0]
    occupied = slots != 0
    is_ribo = np.zeros(slots.shape, dtype=bool)
    is_ribo[occupied] = np.asarray(ribo_IDs, dtype=bool)[slots[occupied].astype(np.int64) - 1]
    has_ribo = is_ribo.any(axis=0)
    needs_move = has_ribo & (membrane[xs, ys, zs] | extracellular[xs, ys, zs])

    if not needs_move.any() and _t127 == 'done':
        # the new centre and cross masks straight from the centre points (same masks as below)
        _pts = np.ascontiguousarray(old_ribo_sites[has_ribo], dtype=np.intp)
        _c, _x = _w127.new_masks(_pts, N_edges[0], N_edges[1], N_edges[2])
        if _RIBO_COMPILED_VERIFY:
            _rc = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool); _rc[xs[has_ribo], ys[has_ribo], zs[has_ribo]] = True
            assert np.array_equal(_rc, _c), 'compiled ribosome pass: centres mismatch'
        ribo_site_dict['ribos']['centers'] = _c
        ribo_site_dict['ribos']['crosses'] = _x
        if len(_pts) == 0:
            ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
        return ribo_site_dict
    if not needs_move.any():
        ribo_center_points[xs[has_ribo], ys[has_ribo], zs[has_ribo]] = True
    else:
        wcm_last_placement_touched_particles = True   # the relocation loop deletes / adds particles
        if plattice is None:
            plattice = lattice.getParticleLatticeView()   # downloads the particles
        for riboSite in old_ribo_sites:
            x, y, z = riboSite[:]
            x_int, y_int, z_int = int(x), int(y), int(z)
        
            if lattice.getOccupancy(z_int, y_int, x_int) > 0:
                parts = getParticlesInSite(plattice, z_int, y_int, x_int)
            
                if len(parts) > 0:
                    for particleIdx in parts:
                        if ribo_IDs[particleIdx-1]:
                            if (membrane[x_int, y_int, z_int] == True) or (extracellular[x_int, y_int, z_int] == True):
                                dist, sheCoordIdx = tree.query([x_int, y_int, z_int])
                                if she_coords is None:
                                    she_coords = _true_coords(cyto_shell)
                                new_position = she_coords[sheCoordIdx]
                            
                                deleteParticle(plattice, z_int, y_int, x_int, particleIdx)
                            
                                new_z, new_y, new_x = int(new_position[2]), int(new_position[1]), int(new_position[0])
                                Occ = lattice.getOccupancy(new_z, new_y, new_x)
                            
                                if Occ < 15:
                                    lattice.addParticle(new_z, new_y, new_x, int(particleIdx))
                                else:
                                    print("Lattice Site Full: ", new_position)
                                    print("Searching for Available Alternative Ribosome Destination Site")
                                
                                    placed = False
                                    place_counter = 1
                                
                                    while not placed:
                                        k_neighbors = min(10 * place_counter, len(she_coords))
                                        siteDists, siteIdxs = tree.query([x_int, y_int, z_int], k=k_neighbors)
                                    
                                        if k_neighbors == 1:
                                            siteIdxs = [siteIdxs]
                                    
                                        for regionCoordIdx in siteIdxs:
                                            test_position = she_coords[regionCoordIdx]
                                            test_z, test_y, test_x = int(test_position[2]), int(test_position[1]), int(test_position[0])
                                        
                                            NewOcc = lattice.getOccupancy(test_z, test_y, test_x)
                                        
                                            if NewOcc < 15:
                                                print("New Destination Site Found: ", test_position)
                                                lattice.addParticle(test_z, test_y, test_x, int(particleIdx))
                                                new_position = test_position
                                                placed = True
                                                break
                                    
                                        place_counter += 1
                                
                                    if not placed:
                                        print("Warning: Could not place ribosome particle", particleIdx)
                                        continue
                            
                                ribo_center_points[int(new_position[0]), int(new_position[1]), int(new_position[2])] = True
                            else:
                                # No relocation needed
                                ribo_center_points[x_int, y_int, z_int] = True
    
    # Reset ribosome site dictionaries
    ribo_site_dict['ribos']['centers'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    ribo_site_dict['ribos']['crosses'] = np.zeros((N_edges[0], N_edges[1], N_edges[2]), dtype=bool)
    
    # Get ribosome center points
    ribo_center_points_coords = _coords(ribo_center_points)
    
    if len(ribo_center_points_coords) == 0:
        return ribo_site_dict
    
    x = ribo_center_points_coords[:, 0]
    y = ribo_center_points_coords[:, 1]
    z = ribo_center_points_coords[:, 2]
    
    ribo_site_dict['ribos']['centers'][x, y, z] = True
    
    x_plus = np.clip(x + 1, 0, N_edges[0] - 1)
    x_minus = np.clip(x - 1, 0, N_edges[0] - 1)
    y_plus = np.clip(y + 1, 0, N_edges[1] - 1)
    y_minus = np.clip(y - 1, 0, N_edges[1] - 1)
    z_plus = np.clip(z + 1, 0, N_edges[2] - 1)
    z_minus = np.clip(z - 1, 0, N_edges[2] - 1)
    
    ribo_site_dict['ribos']['crosses'][x_plus, y, z] = True
    ribo_site_dict['ribos']['crosses'][x_minus, y, z] = True
    ribo_site_dict['ribos']['crosses'][x, y_plus, z] = True
    ribo_site_dict['ribos']['crosses'][x, y_minus, z] = True
    ribo_site_dict['ribos']['crosses'][x, y, z_plus] = True
    ribo_site_dict['ribos']['crosses'][x, y, z_minus] = True

    return ribo_site_dict
#########################################################################################


#########################################################################################
def _put_sites(sv, sites, value, last):
    """Bulk equivalent of ``lattice.setSiteType(site[2], site[1], site[0], value)`` for every row of ``sites``.

    The site-lattice view is indexed [a3, a2, a1] for API arguments (a1, a2, a3), so the call above writes
    ``sv[site[0], site[1], site[2]]``: the view is indexed exactly like the region masks.
    """
    if len(sites):
        sv[sites[:, 0], sites[:, 1], sites[:, 2]] = value
        last[0] = sites[-1]


def updateRiboSites(lattice, ribo_site_dict, region_dict, sim_properties=None):
    """
    Write updated ribosome center/cross masks onto the site lattice.

    ``sim_properties`` is optional (call-site compatibility); unused here.

    Writes go straight into the site-lattice view in the same group order as the former per-site
    ``setSiteType`` loops (later groups overwrite earlier ones, as before). View writes do not clear the
    CUDA lattice's synced flag, so one real ``setSiteType`` call, rewriting a written site's current value,
    is made at the end: the lattice is then uploaded after the hook exactly as before.
    """
    sv = lattice.getSiteLatticeView()
    last = [None]

    # all write groups in one compiled pass (same final site types, same last-written site)
    if len(ribo_site_dict) == 1 and not _RIBO_COMPILED_VERIFY:
        type_dict = next(iter(ribo_site_dict.values()))
        crossID, centerID = type_dict['cross_idx'], type_dict['center_idx']
        R = region_dict
        ms = (R[centerID]['shape'], R[crossID]['shape'], type_dict['centers'], type_dict['crosses'], R['cytoplasm']['shape'],
              R['outer_cytoplasm']['shape'], R['DNA']['shape'], R['extracellular']['shape'], R['membrane']['shape'])
        if _w127_ok(*ms) and sv.dtype == np.uint8 and sv.ndim == 3 and sv.shape == ms[0].shape:
            s = _w127.site_writes(sv, *ms, int(R['cytoplasm']['index']), int(R['outer_cytoplasm']['index']),
                                  int(R['DNA']['index']), int(R['extracellular']['index']), int(R[crossID]['index']), int(R[centerID]['index']))
            R[centerID]['shape'] = type_dict['centers']
            R[crossID]['shape'] = type_dict['crosses']
            if s is not None:
                lattice.setSiteType(int(s[2]), int(s[1]), int(s[0]), int(sv[s[0], s[1], s[2]]))
            return region_dict

    for ribo_type, type_dict in ribo_site_dict.items():
        crossID = type_dict['cross_idx']
        centerID = type_dict['center_idx']

        old_centers = _true_coords(region_dict[centerID]['shape'])

        if len(old_centers) > 0:
            cyto_mask = region_dict['cytoplasm']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]
            outer_cyto_mask = region_dict['outer_cytoplasm']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]
            dna_mask = region_dict['DNA']['shape'][old_centers[:, 0], old_centers[:, 1], old_centers[:, 2]]

            # Priority: cytoplasm > outer_cytoplasm > DNA
            _put_sites(sv, old_centers[cyto_mask], region_dict['cytoplasm']['index'], last)
            _put_sites(sv, old_centers[outer_cyto_mask & ~cyto_mask], region_dict['outer_cytoplasm']['index'], last)
            _put_sites(sv, old_centers[dna_mask & ~cyto_mask & ~outer_cyto_mask], region_dict['DNA']['index'], last)

        old_cross = _true_coords(region_dict[crossID]['shape'])

        if len(old_cross) > 0:
            cyto_mask = region_dict['cytoplasm']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            outer_cyto_mask = region_dict['outer_cytoplasm']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            dna_mask = region_dict['DNA']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]
            extra_mask = region_dict['extracellular']['shape'][old_cross[:, 0], old_cross[:, 1], old_cross[:, 2]]

            # Priority: cytoplasm > outer_cytoplasm > DNA > extracellular
            _put_sites(sv, old_cross[cyto_mask], region_dict['cytoplasm']['index'], last)
            _put_sites(sv, old_cross[outer_cyto_mask & ~cyto_mask], region_dict['outer_cytoplasm']['index'], last)
            _put_sites(sv, old_cross[dna_mask & ~cyto_mask & ~outer_cyto_mask], region_dict['DNA']['index'], last)
            _put_sites(sv, old_cross[extra_mask & ~cyto_mask & ~outer_cyto_mask & ~dna_mask], region_dict['extracellular']['index'], last)

        region_dict[centerID]['shape'] = type_dict['centers']
        region_dict[crossID]['shape'] = type_dict['crosses']

    for ribo_type, type_dict in ribo_site_dict.items():
        crossID = type_dict['cross_idx']
        ribo_sites = _true_coords(type_dict['crosses'])

        if len(ribo_sites) > 0:
            membrane_mask = region_dict['membrane']['shape'][ribo_sites[:, 0], ribo_sites[:, 1], ribo_sites[:, 2]]
            dna_mask = region_dict['DNA']['shape'][ribo_sites[:, 0], ribo_sites[:, 1], ribo_sites[:, 2]]
            _put_sites(sv, ribo_sites[~membrane_mask & ~dna_mask], region_dict[crossID]['index'], last)

    # Set center sites
    for ribo_type, type_dict in ribo_site_dict.items():
        centerID = type_dict['center_idx']
        _put_sites(sv, _true_coords(type_dict['centers']), region_dict[centerID]['index'], last)

    # One real setter call so the CUDA lattice re-uploads (view writes leave its synced flag set).
    if last[0] is not None:
        s = last[0]
        lattice.setSiteType(int(s[2]), int(s[1]), int(s[0]), int(sv[s[0], s[1], s[2]]))

    return region_dict
#########################################################################################


# device-side ribosome lattice reads -------------------------------------------------------------------------------------------------------------------
import os as _wcm_os
_RIBO_DEVICE_VERIFY = _wcm_os.environ.get('WCM_RIBO_DEVICE_VERIFY') is not None


def _wcm_device_reads(lattice):
    """True when the host particles are stale (the solver deferred their download) and the lattice has the device reads."""
    f = getattr(lattice, 'wcmHostParticlesStale', None)
    return bool(f()) if f is not None else False


def _wcm_device_ribo_mask(lattice, ridx, shape):
    """_ribo_sites(plattice, ridx) computed on the device: sites whose packed slots hold particle type ridx."""
    m = np.empty(int(np.prod(shape)), dtype=np.uint8)
    lattice.wcmRiboSiteMask(m, int(ridx))
    return m.view(np.bool_).reshape(shape)


def _wcm_device_slots(lattice, xs, ys, zs, shape):
    """plattice[:, xs, ys, zs, 0] gathered on the device (words x sites, uint32)."""
    flat = np.ravel_multi_index((xs, ys, zs), shape).astype(np.int32)
    words = int(lattice.getMaxOccupancy())
    out = np.empty(words * flat.size, dtype=np.uint32)
    lattice.wcmGatherSlots(flat, out)
    return out.reshape(words, flat.size)


def _wcm_116_check(lattice, ridx, mask, slots, xs, ys, zs):
    """VERIFY: the device mask and slots against the host computation (this downloads the particles)."""
    global _wcm_116_n, _wcm_116_bad
    pl = lattice.getParticleLatticeView()
    ref_mask = _ribo_sites(pl, ridx)
    ref_slots = pl[:, xs, ys, zs, 0]
    same = (ref_mask.shape == mask.shape and bool(np.array_equal(ref_mask, mask)) and ref_slots.dtype == slots.dtype
            and bool(np.array_equal(ref_slots, slots)))
    _wcm_116_n = globals().get('_wcm_116_n', 0) + 1
    _wcm_116_bad = globals().get('_wcm_116_bad', 0) + (0 if same else 1)
    if not same or _wcm_116_n % 1000 == 1:
        print('WCM_RIBO_DEVICE_VERIFY: hook {}: device mask and slots {} the host computation ({} of {} differ)'.format(
            _wcm_116_n, 'IDENTICAL to' if same else 'DIFFERENT from', _wcm_116_bad, _wcm_116_n))
