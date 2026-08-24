"""
Authors: Zane Thornburg

General lattice functions
"""

import numpy as np

# From Tyler Biopolymers
# def deleteParticle(particles, x, y, z, pid):
#     ps = np.array([p for p_ in particles[:,z,y,x,:] for p in p_ if p != pid])
#     pps = 16 # Particles Per Site
# #     pps = cfg.particlesPerSite # Don't know how Tyler added pps as argument to his config
#     ps.resize(pps)
#     ps = ps.reshape((particles.shape[0], particles.shape[4]))
#     particles[:,z,y,x,:] = ps
    
#     return None


from math import floor
from math import log10
def round_sig(x, sig=2):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

    negative = False
    if x < 0:
        negative = True
    x = abs(x)
    if negative:
        return -1*round(x, sig-int(floor(log10(abs(x))))-1)
    elif x==0.0:
        return 0.0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)


def deleteParticle(particles, x, y, z, pid):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

    ps = [p for p_ in particles[:,z,y,x,:] for p in p_ if p != 0]
    if pid in ps:
        ps.remove(pid)
        ps = np.array(ps)
        pps = 16 # Particles Per Site
        ps.resize(pps)
        ps = ps.reshape((particles.shape[0], particles.shape[4]))
        particles[:,z,y,x,:] = ps
    
    return None


def checkParticle(particles, x, y, z, pid):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    """

#     print(x,y,z)
    ps = np.array([p for p_ in particles[:,z,y,x,:] for p in p_ if p == pid])
    if pid in ps:
#         print(ps)
        return True
    else:
        return False
    
    
def getParticlesInSite(particles, x, y, z):
    """
    Inputs:
    Returns:
    Called by:
    Description:
    Optimized version that avoids unnecessary copies while maintaining performance
    for both small and large particle counts per site.
    """
    
    # Get the site slice - shape is (N, 16) where N is first dimension
    site_slice = particles[:, z, y, x, :]
    
    # Use reshape instead of ravel to avoid potential copy
    # Reshape to 1D without copying (if possible) - more efficient than ravel()
    flat_particles = site_slice.reshape(-1)
    
    # Filter non-zero particles - boolean indexing is fast for NumPy arrays
    # This is faster than list comprehension when called many times (like in ribosome code)
    ps = flat_particles[flat_particles != 0]
    
    return ps