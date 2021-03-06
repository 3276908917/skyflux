"""
Proof of concept demonstration. Maybe if we fit the entire chain inside
of one function, we can get some optimizations. First, though,
I want to fix the interpolator. Maybe that will buy us some time,
since we are not iterating through a list of frequencies each time.
"""

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

"""
    Yeah, this is not working. It ices the computer's memory!
    So, instead, I am guessing we do for each source:
        make a 1-source A_tensor
        make a 1-source cold tensor
        (implicitly) garbage collect old 1-source A_tensor
            by relabeling A_tensor to 1-source A_tensor of the next
            source in the sequence
        make a 1-source cold tensor with this new A_tensor,
            then ADD the result to the existing cold tensor
        (implicitly) garbage collect separate 1-source cold tensor
"""
def A_tensor(ra, dec):
    """
    Returned format: an |nu_axis| * |t_axis| * 4 * 4 matrix
    Contains every possible exact A matrix. When performing calculations
    with source index s, time index t, and frequency index f, we say
        this_A = A_tensor[f][s][t]
    """
    global nu_axis
    global t_axis
    
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

"""
To-do:
    1. Saving routine: let us say that it saves our work for every
        frequency completed. That means that, if we crash,
        we lose only 17.21 minutes of work.
        > Laziness bonus: we do not need to work on recovery mechanisms
            (e.g. getting this function to *start* at an arbitrary point)
        Since the full tensor costs about 1.3 MB,
            it will be reasonable to simply save the whole tensor
            (rather than tracking exclusively new work),
            strict upper bound of 336 MB worth of saves.

    two pairs of antannae
        we will use: we want them to be in the same direction, but
        one pair is 14m long and the other is 30m long.

        Answer: 84->85 East to West
                84->86 East to West
            You should calculate the magnitudes of the distances involved.
"""
def cold_tensor(label, ant1, ant2,
                start_index=0, end_index=3871, save_interval=4):
    """
    Returns a giant block of visibility sums. Specifications:
        cold patch: 0 to 8 hours LST in 30 second increments
        full frequency range: 50 to 250 MHz in 1 MHz increments.
    The first return value is the x-axis, also known as the first index.
        It describes the frequency used for that row.
    The second return value is the y-axis, also known as the second index.
        It describes the time used for that column.
    The third return value is the z-axis, also known as the data block.
        It describes the summed visibilities of all ~3000 catalog objects
        for a given time and frequency.
    """
    global nu_axis
    global t_axis
    print("Unix time upon function call:", str(time.time()))
    
    nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
    t_axis = np.arange(0, 2 * np.pi / 3 + MACRO_EPSILON, np.pi / 1440)

    v_tensor = np.zeros((len(nu_axis), len(t_axis), 4), dtype=np.complex128)

    cleaned = demo.cleaned_list()

    percent_interval = 100 / (end_index - start_index + 1)
    percent = 0

    unsaved_counter = 0
    i = start_index
    
    while i < end_index and i < len(cleaned):    
        next_vt = []
        source = cleaned[i]
        
        raI = np.radians(source.ra_angle)
        decI = np.radians(source.dec_angle)
        AI = A_tensor(raI, decI)

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


### Want to write a block that iterates over all possible baselines;
    # ignore duplicates?
    # Create power spectrum and wedge;
        # we don't need to implement the wedge completely;
        # we have HERA pSpec on GitHub
# want a loop over all baselines, over all frequencies,
    # we can use arbitrary fixed LST (within view), and a single source.
# for the wedge, we will want to allow LST to vary over the whole day

def baseline_block(source, LST=None):
    """
    LST none
    """

def load_saves(filename):
    """
    Return a dictionary containing the arrays saved to
    the .npz file at @filename
    """
    a = np.load(filename, allow_pickle=True)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))

# useful for plotting, since I cannot simply pass in op=np.log(np.abs)
logAbs = lambda x: np.log(np.abs(x))

def show(v, op=np.abs):
    plt.imshow(op(v[:, :, 0].T), extent=[50, 250, 8 * 60, 0])
    plt.xlabel('Beam Frequency [MHz]')
    plt.ylabel('LST [minutes]')
    plt.show()
