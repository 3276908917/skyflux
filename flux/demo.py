import matplotlib.pyplot as plt
import numpy as np
import math

from flux import parse

"""
Maybe do a Jupyter notebook with the different histograms
"""

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

def sources_range(start=3, end=5, frq=151):
    assert start < end, "Requested range must be of positive width"
    valid_sources = []
    for gleam_obj in parse.obj_catalog:
        if gleam_obj.flux_by_frq[frq] <= end and \
           gleam_obj.flux_by_frq[frq] >= start:
            valid_sources.append(gleam_obj)
    print("Number of valid sources encountered:", len(valid_sources))
    return valid_sources

# Might be interesting to see whether this changes with different frequencies
def brightest_source(frq=151):
    """
    Return the source with the highest value for integrated flux
    at frequency @frq (in MHz). The source must also satisfy the
    query constaints described in resources/GLEAM_guide.txt.

    There is no error checking to make sure @frq is a valid frequency.
    """
    max_obj = parse.obj_catalog[0]
    for gleam_obj in parse.obj_catalog:
        if gleam_obj.flux_by_frq[frq] > max_obj.flux_by_frq[frq]:
            max_obj = gleam_obj
    print("Largest flux value encountered:", max_obj.flux_by_frq[frq])
    print("Name of associated object:", max_obj.name)
    return max_obj

def is_constrained(value, min_acceptable=None, max_acceptable=None):
    if min_acceptable is not None and value < min_acceptable:
        return False
    if max_acceptable is not None and value > max_acceptable:
        return False
    return True

def hist_data(list_source, frq=151, ln=False, data_lim=None):
    fluxes = []

    if data_lim is not None:
        min_acceptable = data_lim[0]
    else:
        min_acceptable = None
    if data_lim is not None:
        max_acceptable = data_lim[1]
    else:
        max_acceptable = None
        
    for gleam_obj in list_source:
        I = gleam_obj.flux_by_frq[frq]
        if is_constrained(I, min_acceptable, max_acceptable):
            if ln:
                fluxes.append(np.log(I))
            else:
                fluxes.append(I)
                
    return np.array(fluxes)

def brightness_distr(frq=151, ln=False, data_lim=None, ylim=None):
    """
    Generate a histogram for the brigtnesses of all sources
    for a given frequency @frq.
    Where we restrict ourselves to sources between
        @data_lim = (minFlux, maxFlux)
    and we restrict the view according to
        @ylim = (minY, maxY)
    If @ln is True, we take the natural logarithm of each flux before
        plotting it. I found that this helped bring out the appearance
        of the distribution, which was intensely clustered.

    Note: if you want to specify only one boundary (this works
        for both data_lim and ylim), you can simply set the
        irrelevant boundary to None. For example,
            data_lim = (2, None)
        will give all sources brighter than 2 Jy.
    """
    fluxes = hist_data(parse.obj_catalog, frq, ln, data_lim)
    # Naive application of Sturge's Rule to get number of bins
    K = math.ceil(1 + 3.322 * np.log(len(fluxes)))

    fig, ax = frame()
    ax.hist(fluxes, bins=29)

    if ln:
        plt.xlabel("ln([Jy]) at " + str(frq) + " MHz", fontsize=12)
    else:
        plt.xlabel("Flux [Jy] at " + str(frq) + " MHz", fontsize=12)

    plt.ylabel("Frequency", fontsize=12)

    if ylim is not None:
        if ylim[1] is None:
            plt.ylim(bottom=ylim[0])
        elif ylim[0] is None:
            plt.ylim(top=ylim[1])
        else:
            plt.ylim(ylim[0], ylim[1])

"""
76 MHz
Largest flux value encountered: 85.960663
Name of associated object: GLEAM J052257-362727

(This object holds all records on the range [76, 130])

151 MHz
Largest flux value encountered: 55.942432000000004
Name of associated object: GLEAM J052257-362727

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
    for ID in parse.ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(parse.baseline(ant_ID, ID)))

active_ants = list(parse.ant_pos)
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
                  " and " + str(ID2) + " = " + str(parse.baseline(ID1, ID2)))

def phase_factor(ant1, ant2, r, nu=151e6):
    """
    Calculate the phase factor in the direction @r (l, m)
        (we assume that n is of insignificant magnitude)
    and at the frequency @nu
    between two antennae whose ID #s are @ant1 and @ant2.
    When we calculate the baseline (u, v, w), we
        assume that w is of insignificant magnitude.
    """
    b = parse.baseline(ant1, ant2)[0:2] # kill w
    br = np.dot(b, r)
    return np.exp(-2j * np.pi * nu * br / parse.c)
