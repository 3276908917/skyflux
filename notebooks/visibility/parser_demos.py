import matplotlib.pyplot as plt
import numpy as np

import parser

def frame():
    """
    Set up generic infrastructure for an improved-looking plot.
    We return the figure and the axis on which we want to plot.
    """
    fig = plt.figure(figsize = (6, 3))

    plt.subplots_adjust(left=.15, bottom=.2, right=.95, top=.9)
    ax = fig.add_subplot(111)
    
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    return fig, ax

# GLEAMEGCAT section

# Might be interesting to see whether this changes with different frequencies
def brightest_source(frq=151):
    """
    Return the source with the highest value for integrated flux
    at frequency @frq (in MHz). The source must also satisfy the
    query constaints described in resources/GLEAM_guide.txt.

    There is no error checking to make sure @frq is a valid frequency.
    """
    max_obj = parser.obj_catalog[0]
    for gleam_obj in parser.obj_catalog:
        if gleam_obj.flux_by_frq[frq] > max_obj.flux_by_frq[frq]:
            max_obj = gleam_obj
    print("Largest flux value encountered:", max_obj.flux_by_frq[frq])
    print("Name of associated object:", max_obj.name)
    return max_obj

def brightness_distr(frq=151):
    """
    Generate a histogram for the brigtnesses of all sources
    (satisfying the resources/GLEAM_guide.txt constraints)
    for a given frequency.
    """
    fluxes = np.array([
        np.log(gleam_obj.flux_by_frq[frq]) for gleam_obj in parser.obj_catalog
    ])

    fig, ax = frame()
    ax.hist(fluxes, bins=29)
    plt.xlabel("ln([Jy]) at " + str(frq) + " MHz", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

"""
76 MHz
Largest flux value encountered: 37.735298
Name of associated object: GLEAM J084125-754033

(This object holds all records on the range [76, 130])

151 MHz
Largest flux value encountered: 25.322440999999998
Name of associated object: GLEAM J164417-771532

(This object holds all records on the range [143, 227])
"""

# Antenna section

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
    for ID in parser.ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(parser.baseline(ant_ID, ID)))

def all_baselines():
    """
    Print every baseline, without duplicating.

    To-do: eventually I am going to want to figure out how to propagate
        the particular order (in which I am subtracting coordinates)
        out to the integral that I am taking.
    """
    for i in range(len(parser.active_ants)):
        ID1 = parser.active_ants[i]
        for j in range(i + 1, len(parser.active_ants[i + 1:])):
            ID2 = parser.active_ants[j]
            print("Baseline between antennae " + str(ID1) + \
                  " and " + str(ID2) + " = " + str(parser.baseline(ID1, ID2)))

"""
I am not so sure that I really want an integral just yet.
Look at the infinitesimal: it is a solid angle.
Does that not suggest we are integrating over the whole sky?
Currently, our goal is just to pass one point source through the pipe line.
"""

def phase_factor(ant1, ant2, r, nu=151e6):
    """
    Calculate the phase factor in the direction @r (l, m)
        (we assume that n is of insignificant magnitude)
    and at the frequency @nu
    between two antennae whose ID #s are @ant1 and @ant2.
    When we calculate the baseline (u, v, w), we
        assume that w is of insignificant magnitude.
    """
    b = parser.baseline(ant1, ant2)[0:2] # kill w
    br = np.dot(b, r)
    return np.exp(-2j * np.pi * nu * br / parser.c)
