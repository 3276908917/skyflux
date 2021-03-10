import matplotlib.pyplot as plt
import numpy as np
import pickle

from skyflux import catalog
from skyflux import ant
from skyflux import vis
from skyflux import stokes
from skyflux import rot
from skyflux import demo

# disgustingly hacky
MACRO_EPSILON = 0.001

# constants
HOUR = 2 * np.pi / 24
MINUTE = HOUR / 60
SECOND = MINUTE / 60

# we keep these as global parameters to avoid the potential overhead
# of passing by value

# For helices
nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
t_axis = np.arange(0, 24 * HOUR, 30 * SECOND)

nu_rl = range(len(nu_axis))
t_rl = range(len(t_axis))

def A_tensor(ra, dec):
    """
    Returned format: a |nu_axis| * |t_axis| * 4 * 4 matrix
    Contains every possible exact A matrix.
    When performing calculations
    with time index t and frequency index f, we say
        this_A = A_tensor[f][t]
    """
    A_tensor = []

    azs = []
    alts = []
    
    for t in t_axis:
        #! We should definitely vectorize this
        az, alt = rot.eq_to_topo(ra, dec, lst=t, radians=True)
        alts.append(alt)
        azs.append(az)

    for nu in nu_axis:
        J_source = stokes.create_J(
            az=azs, alt=alts, nu=nu, radians=True
        )
        A_source = np.array(
            [stokes.create_A(J=J) for J in J_source]
        )

        A_tensor.append(np.array(A_source))
        
    return np.array(A_tensor)
                   
def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nSimulation: " + percent_status + "% complete.")

def null_source(obj):
    return obj.alpha != obj.alpha
    
def multi_helix(ant1, ant2, sources=catalog.srcs):
    percent_interval = 100 / len(sources)
    percent = 0
    
    helix = single_helix(ant1, ant2, sources[0])
    
    percent += percent_interval
    tick(percent)
    
    for next_obj in sources[1:]:
        if null_source(next_obj):
            continue
        
        next_helix = single_helix(ant1, ant2, next_obj)
        helix = np.add(helix, next_helix)
        
        percent += percent_interval
        tick(percent)
        
    return helix

def single_helix(ant1, ant2, source):
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    s_axis = []
        
    f_layer = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        
        I = vis.get_I(source, nu)
        s = np.array([complex(I), 0, 0, 0])
        A_n = A_full[ni]

        t_layer = []
        for ti in t_rl:
            t = t_axis[ti]

            A_t = A_n[ti]
            
            r = rot.radec2lm(ra, dec, ra0=t)
            phi = ant.phase_factor(ant1, ant2, r, nu)
            
            next_vista = np.dot(np.dot(A_t, s), phi)
            t_layer.append(next_vista)
  
        f_layer.append(np.array(t_layer))

    return np.array(f_layer)
    
def package(block, ptitle):
    """
    Returns a print-ready dictionary,
        particularly for use with the pps.py routines.
        
    Use pickle_dict to save to disk.
    """
    return {'frequencies' : nu_axis,
            'times' : t_axis,
            'picture' : block,
            'title' : ptitle}

def pickle_dict(dict_, label):
    """
    Pickles @dict_ with a file name based on @label.
    While this is a generic routine and will pickle any dictionary,
        the intent is solely for use in conjunction with the
        package_wedge routine.
    """
    with open(label + '.pickle', 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
def auto_helix(ant1, ant2, sources, label, ptitle):
    pickle_dict(package(
        multi_helix(ant1, ant2, sources), ptitle
    ), label)
    
