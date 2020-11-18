import math

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

from skyflux import catalog
from skyflux import ant
from skyflux import vis

# disgustingly hacky

MACRO_EPSILON = 0.001

from skyflux import stokes
from skyflux import rot
from skyflux import demo

# we keep these as global parameters to avoid the potential overhead
# of passing by value
nu_axis = None
t_axis = None

import time

from skyflux import catalog
from skyflux import ant
from skyflux import visibility

### Hard coding, for speed
source = catalog.obj_catalog[3871]
ra = np.radians(source.ra_angle)
dec = np.radians(source.dec_angle)
lst = np.radians(source.ra_angle)

def A_tensor():
    """
    Returned format: a |nu_axis| * |t_axis| * 4 * 4 matrix
    Contains every possible exact A matrix.
    When performing calculations
    with time index t and frequency index f, we say
        this_A = A_tensor[f][t]
    """
    global nu_axis
    global t_axis
    global ra
    global dec
    
    A_tensor = []

    azs = []
    alts = []
    
    for lst in t_axis:
        #! We should definitely vectorize this
        az, alt = rot.eq_to_topo(ra, dec, lst=lst, radians=True)
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

    Returned format: a |nu_axis| * 4 * 4 matrix
    Contains every possible frequency A matrix. When performing calculations
    frequency index f, we say
        this_A = A_tensor[f]
    """
    global nu_axis
    global ra
    global dec
    global lst
    
    A_tensor = []

    az, alt = rot.eq_to_topo(ra, dec, lst=lst, radians=True)

    for nu in nu_axis:
        J_source = stokes.create_J(az=az, alt=alt, nu=nu, radians=True)
        stokes.create_A(J=J_source)

        A_tensor.append(np.array(A_source))
        
    return np.array(A_tensor)

# Scan over all frequencies, for a single source, over all possible baselines
def picture_tensor(source):
    global nu_axis
    global source
    global ra
    global dec
    global lst

    nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
    A = f_only(ra, dec, lst)
    print("\nFinished building partial A tensor.\n")
    
    r = rot.radec2lm(ra, dec, ra0=lst)
    s_axis = []
    
    for ni in range(len(nu_axis)):
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0])

    print("\nFinished building s-vector vector.\n")

    ants = ant.ant_pos.copy()
    outer_ants = ants.copy()
    
    for outer_ant in outer_ants.keys():
        inner_ants = ants.copy()
        del inner_ants[num]

        for inner_ant in inner_ants.keys():
        
            phi = ant.phase_factor(outer_ant, inner_ant, r, nu)
        
            next_vt = []
            for ni in range(len(nu_axis)):
                s = s_axis[ni]
                A_n = A[ni]

                next_vista = np.dot(np.dot(A_n, s), phi)
                next_vt.append(next_vista)

            inner_pos = inner_ants[inner_ant]
            inner_ants[inner_ant] = [inner_pos, next_vt]

        outer_pos = outer_ants[outer_ant]
        outer_ants[outer_ant] = [outer_pos, inner_ants]

    return outer_ants

def wedge_tensor(source):
    global nu_axis
    global source
    global ra
    global dec
    global lst

    nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
    t_axis = np.arange(0, 2 * np.pi, np.pi / 1440)
    
    A_full = A_tensor()
    print("\nFinished building A tensor.\n")

    r = rot.radec2lm(ra, dec, ra0=lst)
        s_axis = []
    
    for ni in range(len(nu_axis)):
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0])

    print("\nFinished building s-vector vector.\n")

    ants = ant.ant_pos.copy()
    outer_ants = ants.copy()

    percent_interval = 100 / (len(ants) * (len(ants) - 1))
    percent = 0

    for outer_ant in outer_ants.keys():
        inner_ants = ants.copy()
        del inner_ants[num]

        for inner_ant in inner_ants.keys():
            phi = ant.phase_factor(outer_ant, inner_ant, r, nu)
        
            next_vt = []
            for ni in range(len(nu_axis)):
                nu = nu_axis[ni]
                next_vt.append([])

                I = vis.get_I(source, nu)
                s = np.array([complex(I), 0, 0, 0])

                A_n = A_full[ni]

                for ti in range(len(t_axis)):
                    t = t_axis[ti]

                    A = A_n[ti]
                    
                    next_vista = np.dot(np.dot(A, s), phi)#####

                next_vista = np.dot(np.dot(A_n, s), phi)
                next_vt.append(next_vista)

            inner_pos = inner_ants[inner_ant]
            inner_ants[inner_ant] = [inner_pos, next_vt]

        outer_pos = outer_ants[outer_ant]
        outer_ants[outer_ant] = [outer_pos, inner_ants]

    return outer_ants

def cold_tensor(label, ant1, ant2,
                start_index=0, end_index=3871, save_interval=4):

        for ni in range(len(nu_axis)):
            nu = nu_axis[ni]
            next_vt.append([])

            I = vis.get_I(source, nu)
            s = np.array([complex(I), 0, 0, 0])

            A_n = AI[ni]
            
            for ti in range(len(t_axis)):
                t = t_axis[ti]

                A = A_n[ti]
                r = rot.radec2lm(raI, decI, ra0=t)
                phi = ant.phase_factor(ant1, ant2, r, nu)

                next_vista = np.dot(np.dot(A, s), phi)
                next_vt[len(next_vt) - 1].append(next_vista)

        v_tensor += np.array(next_vt)
        
        percent += percent_interval
        percent_status = str(np.around(percent, 4))
        print("Visibility tensor: " + percent_status + \
              "% complete (finished i=" + str(i) + ").")

        unsaved_counter += 1
        if unsaved_counter > save_interval:
            np.savez("backup_" + label, na=nu_axis, ta=t_axis, vt=v_tensor,
                     dying_index=np.array(i))
            unsaved_counter = 0

        i += 1

    np.savez("backup_" + label, na=nu_axis, ta=t_axis, vt=v_tensor,
                     dying_index=np.array(-1))    

"""
Visualization:
plt.imshow(np.abs(vt[:, :, 0].T), extent=[50, 250, 8 * 60, 0])
plt.xlabel('Beam Frequency [MHz]')
plt.ylabel('LST [minutes]')
plt.show()
"""
