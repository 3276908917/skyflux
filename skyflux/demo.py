"""
Utilities for displaying various plots concerning the distributions
of the GLEAM objects in the parser buffer (catalog.obj_catalog).
"""

import matplotlib.pyplot as plt
import numpy as np
import math

from skyflux import catalog
from skyflux import ant

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

# Visibility section

if False: # keep this section on ice until we are satisfied
    nside = 128
    theta, phi = hp.pix2ang(nside, np.arange(12 * nside * nside))
    az = phi
    alt = np.pi/2 - theta

    J = stokes.create_J(az=az, alt=alt, radians=True)
    #! Does this work the same, if I insert az/alt directly into create_A?
    A = np.array([stokes.create_A(J=Ji) for Ji in J])

    def project_J(J, data_transform=np.abs, rep=hp.orthview):
        def put_subplot(i, j, panel, ttl=None):
            if ttl is None:
                ttl = str(i) + ', ' + str(j)
            # hard-coding is always bad
            if rep is hp.cartview:
                rep(data_transform(J[:, i, j]),
                    half_sky=True, sub=[3, 2, panel], title=ttl)
            else:
                rep(data_transform(J[:, i, j]), rot=[0, 90],
                    half_sky=True, sub=[3, 2, panel], title=ttl)

        put_subplot(0, 0, 1, 'xx')
        put_subplot(1, 0, 2, 'yx')
        put_subplot(0, 1, 3, 'xy')
        put_subplot(1, 1, 4, 'yy')

    def project_A(A, data_transform=np.abs, rep=hp.orthview):
        def put_subplot(i, j, panel, ttl=None):
            if ttl is None:
                ttl = str(i) + ', ' + str(j)
            if rep is hp.cartview:
                rep(data_transform(A[:, i, j]),
                    half_sky=True, sub=[5, 4, panel], title=ttl)
            else:
                rep(data_transform(A[:, i, j]), rot=[0, 90]
                    half_sky=True, sub=[5, 4, panel], title=ttl)
        
        put_subplot(0, 0, 1, 'I\' <- I')
        put_subplot(0, 1, 2)
        put_subplot(0, 2, 3)
        put_subplot(0, 3, 4)
        
        put_subplot(1, 0, 5)
        put_subplot(1, 1, 6, 'Q\' <- Q')
        put_subplot(1, 2, 7)
        put_subplot(1, 3, 8)
        
        put_subplot(2, 0, 9)
        put_subplot(2, 1, 10)
        put_subplot(2, 2, 11, 'U\' <- U')
        put_subplot(2, 3, 12)
        
        put_subplot(3, 0, 13)
        put_subplot(3, 1, 14)
        put_subplot(3, 2, 15)
        put_subplot(3, 3, 16, 'V\' <- V')

# GLEAMEGCAT section

def sources_range(start=3, end=5, frq=151):
    """
    Return all sources in the current parse buffer
        (see obj_catalog in catalog.py)
    for which the flux at the frequency @frq (in MHz) is
        at least @start [Jy] and
        at most @end [Jy]
    This function also prints the number of such sources encountered.
    """
    assert start < end, "Requested range must be of positive width"
    valid_sources = []
    for gleam_obj in catalog.obj_catalog:
        if gleam_obj.flux_by_frq[frq] <= end and \
           gleam_obj.flux_by_frq[frq] >= start:
            valid_sources.append(gleam_obj)
    print("Number of valid sources encountered:", len(valid_sources))
    return valid_sources

def brightest_source(frq=151):
    """
    Return the source with the highest value for integrated flux
        at frequency @frq (in MHz).
    The source must also satisfy the
    query constaints described in resources/GLEAM_guide.txt.

    There is no error checking to make sure @frq is a valid frequency.
    """
    max_obj = catalog.obj_catalog[0]
    for gleam_obj in catalog.obj_catalog:
        if gleam_obj.flux_by_frq[frq] > max_obj.flux_by_frq[frq]:
            max_obj = gleam_obj
    print("Largest flux value encountered:", max_obj.flux_by_frq[frq])
    print("Name of associated object:", max_obj.name)
    return max_obj

def is_constrained(value, min_acceptable=None, max_acceptable=None):
    """
    Return True if @max_acceptable >= @value >= @min_acceptable,
        False otherwise.
    Although this is written as a general helper function,
        it is chiefly used in filtering the data for the histogram.
    """
    if min_acceptable is not None and value < min_acceptable:
        return False
    if max_acceptable is not None and value > max_acceptable:
        return False
    return True

def hist_data(list_source, frq=151, ln=False, data_lim=None):
    """
    For every GLEAM_entry object
        as defined in catalog.py
    in @list_source,
    extract its flux at frequency @frq (in MHz)

    Return an array containing those fluxes which
    satisfy the boundaries established by
        @data_lim = (@min_acceptable_value, @max_acceptable_value).
    To leave an extreme unbounded, enter None as that value.
    
    If @ln = True, we also take the natural logarithm of all the
    fluxes that we return. This is helpful, although not rigorous,
    for observing some of the more exponential spread patterns,
    such as the radio flux distribution as 151 MHz.
    """
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
    for a given frequency @frq (in MHz).
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
    fluxes = hist_data(catalog.obj_catalog, frq, ln, data_lim)
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
    for ID in ant.ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(ant.baseline(ant_ID, ID)))

active_ants = list(ant.ant_pos)
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
                  " and " + str(ID2) + " = " + str(ant.baseline(ID1, ID2)))
