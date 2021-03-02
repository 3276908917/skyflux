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
hour = 2 * np.pi / 24
minute = hour / 60

### Hard coding, for speed ###
source = catalog.obj_catalog[3871]
lst0 = np.radians(source.ra_angle)

# we keep these as global parameters to avoid the potential overhead
# of passing by value

# For wedges
# nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 4e6)
# For helices
nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)

nu_rl = range(len(nu_axis))

# For wedges
# t_axis = np.arange(lst0 - hour, lst0 + hour, 4 * minute)
# For helices: 30 second intervals
t_axis = np.arange(0, 2 * np.pi, np.pi / 1440)
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
        J_source = stokes.create_J(az=azs, alt=alts, nu=nu, radians=True)
        A_source = np.array([stokes.create_A(J=J) for J in J_source])

        A_tensor.append(np.array(A_source))
        
    return np.array(A_tensor)
                   
def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nSimulation: " + percent_status + "% complete.")

def null_source(obj):
    return obj.alpha != obj.alpha

def single_helix(ant1, ant2, source):
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    r = rot.radec2lm(ra, dec, ra0=lst0)
    s_axis = []
        
    f_layer = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        phi = ant.phase_factor(ant1, ant2, r, nu)
        
        I = vis.get_I(source, nu)
        s = np.array([complex(I), 0, 0, 0])
        A_n = A_full[ni]

        t_layer = []
        for ti in t_rl:
            t = t_axis[ti]

            A_t = A_n[ti]
            
            next_vista = np.dot(np.dot(A_t, s), phi)
            t_layer.append(next_vista)
  
        f_layer.append(np.array(t_layer))

    return np.array(f_layer)
    
