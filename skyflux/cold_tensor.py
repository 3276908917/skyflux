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

def A_tensor(nu_axis, t_axis, sources):
    """
    Returned format: an |nu_axis| * |sources| * |t_axis| * 4 * 4 matrix
    Contains every possible exact A matrix. When performing calculations
    with source index s, time index t, and frequency index f, we say
        this_A = A_tensor[f][s][t]
    """
    A_tensor = []
    for nu in nu_axis:
        A_tensor.append([])
        for source in sources:
            ra = np.radians(sources.ra_angle)
            dec = np.radians(souces.dec_angle)
            azs = []
            alts = []
            for lst in t_axis:
                az, alt = eq_to_topo(ra, dec, lst=lst, radians=True)
                alts.append(alt)
                azs.append(az)
            J_source = stokes.create_J(az=azs, alt=alts, nu=nu, radians=True)
            A_source = np.array([stokes.create_A(J) for J in J_source])

            A_tensor[len(A_tensor) - 1].append(A_source)
    return A_tensor

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
    2. Extrapolation routine:
        We want to fix the code, which for some reason uses a linear
            interpolation *on top* of a power law relation.
        The fix should be simple, but we will certainly want to re-run
        the single_source_over_nu notebook to make sure our results
        did not change.
    3. Examine the antennae chart and decide which two pairs of antannae
        we will use: we want them to be in the same direction, but
        one pair is 14m long and the other is 30m long.
"""
def cold_tensor(ant1, ant2, sources=None):
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
    import time as t
    print("Unix time upon function call:", str(t.time()))

    percent_interval = 100 / 961 / 251
    
    nu_axis = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)
    t_axis = np.arange(0, 2 * np.pi / 3 + MACRO_EPSILON, np.pi / 1440)
    v_tensor = []

    if sources is None:
        sources = catalog.obj_catalog.copy()


    ### We want to pre-generate all the A matrices

    A_tensor = []

    for nu in nu_axis

    ###

    percent = 0
    for nu in nu_axis:
        v_tensor.append([])
        for t in t_axis:
            next_vista = np.array([0j, 0j, 0j, 0j])
            for source in sources:
                next_vista += vis.visibility(ant1, ant2, source, nu=nu, time=t)

            percent += percent_interval
            v_tensor[len(v_tensor) - 1].append(next_vista)

            percent_status = str(np.around(percent, 4))
            print("Visibility tensor: " + percent_status + "% complete.")

    return nu_axis, t_axis, np.array(v_tensor)
