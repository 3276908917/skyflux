"""
hera_pspec:
    spectral window can be [(145 MHz, 155 MHz)]

    First channel is 50 MHz which corresponds to index 0
        like, ch = np.arange(201)
    indices into frequency array
"""

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

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
nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 4e6)
nu_rl = range(len(nu_axis))

t_axis = np.arange(lst0 - hour, lst0 + hour, 4 * minute)
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

def f_only():
    """
    I am only accepting lst in radians.

    Counterpart to above function, but used when only a single
        LST is of concern.

    Returned format: a |nu_axis| * 4 * 4 matrix
    Contains every possible frequency A matrix. When performing calculations
    frequency index f, we say
        this_A = A_tensor[f]
    """
    A_tensor = []

    az, alt = rot.eq_to_topo(ra, dec, lst=lst0, radians=True)

    for nu in nu_axis:
        J_source = stokes.create_J(az=az, alt=alt, nu=nu, radians=True)
        A_source = stokes.create_A(J=J_source)

        A_tensor.append(np.array(A_source))
        
    return np.array(A_tensor)

import pickle

def pickle_dict(dict_, label):
    with open(label + '.pickle', 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)

# super_dict = picture_tensor()

# Scan over all frequencies, for a single source, over all possible baselines
def picture_tensor():
    raise NotImplementedError("Still processes just one source.")

    nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
    A = f_only()
    print("\nFinished building partial A-tensor.\n")
    
    r = rot.radec2lm(ra, dec, ra0=lst0)
    s_axis = []
    
    for ni in range(len(nu_axis)):
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0]))

    print("\nFinished building s-vector vector.\n")

    ants = ant.ant_pos.copy()
    outer_ants = ants.copy()
    
    for outer_ant in outer_ants.keys():
        inner_ants = ants.copy()
        del inner_ants[outer_ant]

        for inner_ant in inner_ants.keys():
        
            phi = ant.phase_factor(outer_ant, inner_ant, r, nu)
        
            next_vt = []
            for ni in range(len(nu_axis)):
                s = s_axis[ni]
                A_n = A[ni]

                next_vista = np.dot(np.dot(A_n, s), phi)
                next_vt.append(next_vista)

            inner_ants[inner_ant] = next_vt

        outer_ants[outer_ant] = inner_ants

    return outer_ants

def merge_wedges(wedge1, wedge2):
    """ We assume that both wedges have the same format:
        ant1, ant2, nu, t hierarchies are exactly the same.
        
        WARNING: this is not write-safe!
        To conserve memory, we modify the wedge2 parameter."""

    for ant1 in wedge1.keys():
        for ant2 in wedge1[ant1].keys():
            for nu_idx in nu_rl:
                system = wedge1[ant1][ant2]
                for t_idx in t_rl:
                    visibility2 = wedge2[ant1][ant2][nu_idx][t_idx]
                    system[t_idx] += visibility2
                   

def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nWedge tensor: " + percent_status + "% complete.")

def full_wedge():
    
    percent_interval = 100 / len(catalog.obj_catalog)
    percent = 0
    
    wedge = single_wedge(catalog.obj_catalog[0])
    
    percent += percent_interval
    tick(percent)
    
    for next_obj in catalog.obj_catalog[1:]:
        next_wedge = single_wedge(next_obj)
        merge_wedges(wedge, next_wedge)

        percent += percent_interval
        tick(percent)
        
    return wedge

outer_ants = ant.ant_pos.copy()

def single_wedge(source):
    # note
    # it does not really make sense to do a two-hour interval
    # with a full sky,
    # but I do not want to risk large file sizes,
    # so the simulation will be deliberately truncated...
    
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    r = rot.radec2lm(ra, dec, ra0=lst0)
    s_axis = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0]))
   
    for outer_ant in outer_ants.keys():
        inner_ants = outer_ants.copy()
        del inner_ants[outer_ant]

        for inner_ant in inner_ants.keys():
        
            phi = ant.phase_factor(outer_ant, inner_ant, r, nu)
        
            f_layer = []
            for ni in nu_rl:
                s = s_axis[ni]
                A_n = A_full[ni]

                t_layer = []
                for ti in t_rl:
                    t = t_axis[ti]

                    A_t = A_n[ti]
                    
                    next_vista = np.dot(np.dot(A_t, s), phi)
                    t_layer.append(next_vista)
          
                f_layer.append(t_layer)

            inner_ants[inner_ant] = f_layer

        outer_ants[outer_ant] = inner_ants

    return outer_ants

