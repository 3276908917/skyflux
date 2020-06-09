import pickle
import glob
import numpy as np
import time
import astropy
import astropy.time

import rotations

# Antenna section

ant_pos = dict(pickle.load(open("ant_dict.pk", "rb")))

"""
b = (u, v, w) is the
vector representing the coordinates in meters in the plane
of the array
"""

def list_baselines(ant_ID):
    """
    Print every baseline that features antenna # @ant_ID
    """
    print ("Baselines between antenna " + str(ant_ID) + " and antenna...")
    for ID in ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(baseline(ant_ID, ID)))

active_ants = list(ant_pos)
active_ants.sort()

def all_baselines():
    """
    Print every baseline, without duplicating.

    To-do: eventually I am going to want to figure out how to propagate
        the particular order (in which I am subtracting coordinates)
        out to the integral that I am taking.
    """
    for i in range(len(active_ants)):
        ID1 = active_ants[i]
        for j in range(i + 1, len(active_ants[i + 1:])):
            ID2 = active_ants[j]
            print("Baseline between antennae " + str(ID1) + \
                  " and " + str(ID2) + " = " + str(baseline(ID1, ID2)))

"""
I am not so sure that I really want an integral just yet.
Look at the infinitesimal: it is a solid angle.
Does that not suggest we are integrating over the whole sky?
Currently, our goal is just to pass one point source through the pipe line.
"""

S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

c = 299792458 # m / s

def phase_factor(ant1, ant2, r, nu=151e6):
    """
    Calculate the phase factor in the direction @r (l, m)
        (we assume that n is of insignificant magnitude)
    and at the frequency @nu
    between two antennae whose ID #s are @ant1 and @ant2.
    When we calculate the baseline (u, v, w), we
        assume that w is of insignificant magnitude.
    """
    b = baseline(ant1, ant2)[0:2] # kill w
    br = np.dot(b, r)
    return np.exp(-2j * np.pi * nu * br / c)

def visibility_integrand(J, source, nu=151e6):
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    r = rotations.raddec2lm(ra, dec)
    
    J_outer = np.kron(J, np.conj(J))
    A = np.dot(np.dot(np.linalg.inv(S), J_outer), S)

    # Do I want to sum for all possible baselines?
    # If not: should I make a function to evaluate the integrand for a single baseline?
    phi = phase_sum(r, nu)
    return np.dot(np.dot(A, s), phi)

"""
Sum your sources, not your baselines.
Summation over different values r-hat

Integration is a summation over sources.

(I think the output should be a scalar, not a 4x1)

The graph of times and positions

A will change according to time. Equation 3 assumes that we have
just one integration time.

The change in A over time represents a drift scan,
like the Gaussians that you recently visualized.

Time you can do, you have all the information

use constant flux independent of frequency
(this is the same as saying that alpha is equal to one)
then varying frequency from 100-200 MHz
'assume it is the same beam between 100 and 200 MHz'
    beam variation is what we will eventually be handling anyway

15 m or 30 m baseline East-West

I need to plot the A matrix, to see the leakage terms.
Plot it in healpix.

(Ask for her power spectrum notebook plot for a single baseline,
if you complete all of your work.)
"""

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
