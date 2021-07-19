"""
Utilities for displaying various plots concerning the distributions
of the GLEAM objects in the parser buffer (catalog.obj_catalog).
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

# general helpers and utilities

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

# Visibility section

"""
Some ideas for testing the output of vis_tensor

src = sf.catalog.obj_catalog[0]
single_nu, single_t, v_tensor = sf.demo.vis_tensor(23, 37, np.array([src]))

# remember that each element of the visibility tensor is
# a 4x1 complex vector of I, Q, U, and V
plt.plot(single_nu, v_tensor[:, 0, 0])
plt.plot(single_t, v_tensor[0, :, 0])

# begin actual experiment
nu14, lst14, vt14 = vis_tensor(23, 37)
nu30, lst30, vt30 = vis_tensor(37, 68)

ss = sf.catalog.obj_catalog[3871]
nu14_s, lst14_s, vt14_s = vis_tensor(23, 37, ss)
nu30_s, lst30_s, vt30_s = vis_tensor(37, 68, ss)

np.savez_compressed('14m_baseline',
    frequency_axis=nu14, lst_axis=lst14, visibility_tensor=vt14)
np.savez_compressed('30m_baseline',
    frequency_axis=nu30, lst_axis=lst30, visibility_tensor=vt30)

np.savez_compressed('14m_single_source',
    frequency_axis=nu14_s, lst_axis=lst14_s, visibility_tensor=vt14_s)
np.savez_compressed('30m_single_source',
    frequency_axis=nu30_s, lst_axis=lst30_s, visibility_tensor=vt30_s)
"""
def vis_tensor(ant1, ant2, sources=None):
    """
    Returns a giant block of visibility sums.
    The first return value is the x-axis, also known as the first index.
        It describes the frequency used for that row.
    The second return value is the y-axis, also known as the second index.
        It describes the time used for that column.
    The third return value is the z-axis, also known as the data block.
        It describes the summed visibilities of all ~3000 catalog objects
        for a given time and frequency.

    def vis_tensor(ant1, ant2, sources=None):

        import time as t
        print("Unix time upon function call:", str(t.time()))

        percent_interval = 100 / 144 / 150
    
        nu_axis = np.arange(77e6, 226e6 + sf.demo.MACRO_EPSILON, 1e6)
        t_axis = np.arange(0, 2 * np.pi, np.pi / 72)
        v_tensor = []

        if sources is None:
            sources = sf.catalog.obj_catalog.copy()

        percent = 0
        for nu in nu_axis:
            v_tensor.append([])
            for t in t_axis:
                next_vista = np.array([0j, 0j, 0j, 0j])
                for source in sources:
                    next_vista += sf.vis.visibility(
                    ant1, ant2, source, nu=nu, time=t)

                percent += percent_interval
                v_tensor[len(v_tensor) - 1].append(next_vista)

                percent_status = str(np.around(percent, 4))
                print("Visibility tensor: " + percent_status + "% complete.")

        return nu_axis, t_axis, np.array(v_tensor)
    """
    import time as t
    print("Unix time upon function call:", str(t.time()))

    percent_interval = 100 / 144 / 150
    
    nu_axis = np.arange(77e6, 226e6 + MACRO_EPSILON, 1e6)
    t_axis = np.arange(0, 2 * np.pi, np.pi / 72)
    v_tensor = []

    if sources is None:
        sources = catalog.obj_catalog.copy()

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

### todo: we want a command that will force all of the scales to run from the same values

# Appearances: jones_matrices/A_Catalog.ipynb
def project_J(J, data_transform=np.abs, rep=hp.orthview,
              widest_scale=False, max_=None, min_=None):
    """
    Generate a 2x2 plot of the four Jones components,
    assuming that J has conventional formatting [[xx, xy], [yx, yy]]
    @data_transform : function applied to each index of the J matrix
    @rep : one of the healpy projection functions, such as
        orthview, cartview, and mollview
    """
    scale_max = max_
    scale_min = min_

    if widest_scale:
        if max_ is not None or min_ is not None:
            raise TypeError("Cannot have automatic and manual " + \
                            "axes simultaneously.")
        scale_max = np.max(J)
        scale_min = np.min(J)
    
    def put_subplot(i, j, panel, ttl=None):
        if ttl is None:
            ttl = str(i) + ', ' + str(j)
        # hard-coding is always bad
        if rep is hp.cartview:
            rep(data_transform(J[:, i, j]),
                half_sky=True, sub=[3, 2, panel], title=ttl,
                max=scale_max, min=scale_min)
        else:
            rep(data_transform(J[:, i, j]), rot=[0, 90],
                half_sky=True, sub=[3, 2, panel], title=ttl,
                max=scale_max, min=scale_min)

    put_subplot(0, 0, 1, 'xx')
    put_subplot(1, 0, 2, 'yx')
    put_subplot(0, 1, 3, 'xy')
    put_subplot(1, 1, 4, 'yy')

# Appearances: jones_matrices/A_Catalog.ipynb
def project_A(A, data_transform=np.abs, rep=hp.orthview,
              widest_scale=False, max_=None, min_=None):
    """
    Generate a 4x4 plot of the sixteen Mueller components,
    assuming that A has conventional formatting
        [[I' <- I, I' <- Q, I' <- U, I' <- V],
        [Q' <- I, Q' <- Q, Q' <- U, Q' <- V],
        [U' <- I, U' <- Q, U' <- U, U' <- V],
        [V' <- I, V' <- Q, V' <- U, V' <- V]]
    @data_transform : function applied to each index of the J matrix
    @rep : one of the healpy projection functions, such as
        orthview, cartview, and mollview
    """
    scale_max = max_
    scale_min = min_

    # Actually, I think it is legal to specify a min without a max
    #if scale_max is None and scale_min is not None or \
    #   scale_max is not None and scale_min is None:
    #    raise TypeError("")

    if widest_scale:
        if max_ is not None or min_ is not None:
            raise TypeError("Cannot have automatic and manual " + \
                            "axes simultaneously.")
        scale_max = np.max(A)
        scale_min = np.min(A)
    
    def put_subplot(i, j, panel, ttl=None):
        if ttl is None:
            ttl = str(i) + ', ' + str(j)
        if rep is hp.orthview:
            rep(data_transform(A[:, i, j]), rot=[0, 90],
                half_sky=True, sub=[5, 4, panel], title=ttl,
                max=scale_max, min=scale_min)
        else:
            rep(data_transform(A[:, i, j]), rot=[0, 90],
                sub=[5, 4, panel], title=ttl,
                max=scale_max, min=scale_min)

    put_subplot(0, 0, 1, 'I\' <- I')
    put_subplot(0, 1, 2, 'I\' <- Q')
    put_subplot(0, 2, 3, 'I\' <- U')
    put_subplot(0, 3, 4, 'I\' <- V')

    put_subplot(1, 0, 5, 'Q\' <- I')
    put_subplot(1, 1, 6, 'Q\' <- Q')
    put_subplot(1, 2, 7, 'Q\' <- U')
    put_subplot(1, 3, 8, 'Q\' <- V')

    put_subplot(2, 0, 9, 'U\' <- I')
    put_subplot(2, 1, 10, 'U\' <- Q')
    put_subplot(2, 2, 11, 'U\' <- U')
    put_subplot(2, 3, 12, 'U\' <- V')

    put_subplot(3, 0, 13, 'V\' <- I')
    put_subplot(3, 1, 14, 'V\' <- Q')
    put_subplot(3, 2, 15, 'V\' <- U')
    put_subplot(3, 3, 16, 'V\' <- V')

# GLEAMEGCAT section

def lookup(name):
    """
    Given the @name field for a GLEAMEGCAT object,
    return the index (in the master array) corresponding to that object.

    Returns -1 if the name was not found in the catalog.
    """
    for i in range(len(catalog.obj_catalog)):
        if catalog.obj_catalog[i].name == name:
            return i
    return -1

def cleaned_list():
    """
    Returns a copy of the GLEAMEGCAT list,
    for which we have eliminated all entries with
    unspecified spectral indices. By construction,
    objects without spectral indices were
    automatically assigned 'NaN' for the field.
    """
    ws_oc = catalog.srcs.copy() # write-safe read copy for
        # the GLEAM object catalog
    cat = catalog.srcs.copy()
    # we loop in reverse, to avoid concurrent mod. exceptions
    for i in range(len(ws_oc) - 1, 0, -1):
        # classic. The easiest way to check if a value is NaN:
        # it won't equal itself
        if ws_oc[i].alpha != ws_oc[i].alpha:
            cat = np.delete(cat, i)
    return cat

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

def brightest_source(frq=151, sliced_list=catalog.srcs):
    """
    Return the source with the highest value for integrated flux
        at frequency @frq (in MHz).
    The source must also satisfy the
    query constaints described in resources/GLEAM_guide.txt.

    There is no error checking to make sure @frq is a valid frequency.
    """
    max_obj = sliced_list[0]
    for gleam_obj in sliced_list:
        if gleam_obj.flux_by_frq[frq] > max_obj.flux_by_frq[frq]:
            max_obj = gleam_obj
    print("Largest flux value encountered:", max_obj.flux_by_frq[frq])
    print("Name of associated object:", max_obj.name)
    print("Index of associated object:", lookup(max_obj.name))
    return max_obj

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

# Antenna section

"""
b = (u, v, w) is the
vector representing the coordinates in meters in the plane
of the array
"""

def list_baselines(ant_ID):
    """
    Print every available baseline that features antenna # @ant_ID
    """
    print ("Baselines between antenna " + str(ant_ID) + " and antenna...")
    for ID in ant.ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(ant.baseline(ant_ID, ID)))

active_ants = list(ant.ant_pos)
active_ants.sort()

def all_baselines():
    """
    Print every available baseline, without duplicating.
    """
    for i in range(len(active_ants)):
        ID1 = active_ants[i]
        for j in range(i + 1, len(active_ants[i + 1:])):
            ID2 = active_ants[j]
            print("Baseline between antennae " + str(ID1) + \
                  " and " + str(ID2) + " = " + str(ant.baseline(ID1, ID2)))

def find_baseline(objective):
    """
    Returns the baseline
        (antenna 1 ID #, antenna 2 ID #)
    that is closest in length to
    objective out of all possible antenna pairs.
    """
    best_err = float('inf')
    best_ants = None
    
    for i in range(len(active_ants)):
        ID1 = active_ants[i]
        for j in range(i + 1, len(active_ants[i + 1:])):
            ID2 = active_ants[j]
            b = ant.baseline(ID1, ID2)
            sq_err = (objective - np.linalg.norm(b)) ** 2

            if sq_err < best_err:
                best_err = sq_err
                best_ants = (ID1, ID2)

    print("Best squared error:", best_err)
    return best_ants
