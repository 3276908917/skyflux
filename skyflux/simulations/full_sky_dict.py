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
HOUR = 2 * np.pi / 24
MINUTE = HOUR / 60
SECOND = MINUTE / 60

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
t_axis = np.arange(0, 24 * HOUR, 30 * SECOND)
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
    raise NotImplementedError("Still processes just one source.")
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
        To conserve memory, we modify the wedge1 parameter."""


    sum_ = {}
    for ant1 in wedge1.keys():
        sum_[ant1] = {}
        for ant2 in wedge1[ant1].keys():
            sum_[ant1][ant2] = np.zeros(
                (len(nu_rl), len(t_rl), 4), dtype=np.complex128
            )
            for nu_idx in nu_rl:
                for t_idx in t_rl:
                    vis1 = wedge1[ant1][ant2][nu_idx][t_idx]
                    vis2 = wedge2[ant1][ant2][nu_idx][t_idx]
                    
                    """ Debugging block
                    print(vis1)
                    print()
                    print(vis2)
                    print()
                    print(sum_[ant1][ant2][nu_idx][t_idx])
                    """
                    
                    sum_[ant1][ant2][nu_idx][t_idx] = \
                        np.add(vis1, vis2)
                    if np.array_equal(
                        sum_[ant1][ant2][nu_idx][t_idx],
                        np.zeros(4)
                    ):
                        print("\nZero encountered!")
                        print(vis1)
                        print(vis2)
                        print()
                
    return sum_
                
    """
    for ant1 in wedge1.keys():
        for ant2 in wedge1[ant1].keys():
            for nu_idx in nu_rl:
                system = wedge1[ant1][ant2][nu_idx]
                for t_idx in t_rl:
                    visibility2 = wedge2[ant1][ant2][nu_idx][t_idx]
                    #print("\n" + str(type(visibility2)) + "\n")
                    wedge1[ant1][ant2][nu_idx][t_idx] = np.add(
                        visibility2, system[t_idx])
                    if np.array_equal(system[t_idx], np.zeros(4)):
                        print("\nZero encountered!")
                        print(system[t_idx] - visibility2)
                        print(visibility2)
                        print()
    """
                   
def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nSimulation: " + percent_status + "% complete.")

def null_source(obj):
    return obj.alpha != obj.alpha

def full_wedge(list_sources=catalog.srcs):
    percent_interval = 100 / len(list_sources)
    percent = 0
    
    wedge = single_wedge(list_sources[0])
    
    percent += percent_interval
    tick(percent)
    
    for next_obj in list_sources:
        if null_source(next_obj):
            continue
        
        next_wedge = single_wedge(next_obj)
        wedge = merge_wedges(wedge, next_wedge)
        
        # Have we successfully moved beyond this bug?
        #print(wedge[136][140])

        percent += percent_interval
        tick(percent)
        
    return wedge
    
def multi_helix(ant1, ant2, sources=catalog.srcs):
    percent_interval = 100 / len(sources)
    percent = 0
    
    helix = single_helix(ant1, ant2, sources[0])
    
    percent += percent_interval
    tick(percent)
    
    for next_obj in sources[1:]:
        if null_source(next_obj):
            continue
        
        next_helix = single_helix(ant1, ant2, next_obj)
        helix = np.add(helix, next_helix)
        
        percent += percent_interval
        tick(percent)
        
    return helix

outer_ants = ant.ant_pos.copy()

def single_helix(ant1, ant2, source):
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    s_axis = []
        
    f_layer = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        
        I = vis.get_I(source, nu)
        s = np.array([complex(I), 0, 0, 0])
        A_n = A_full[ni]

        t_layer = []
        for ti in t_rl:
            t = t_axis[ti]

            A_t = A_n[ti]
            
            r = rot.radec2lm(ra, dec, ra0=t)
            phi = ant.phase_factor(ant1, ant2, r, nu)
            
            next_vista = np.dot(np.dot(A_t, s), phi)
            t_layer.append(next_vista)
  
        f_layer.append(np.array(t_layer))

    return np.array(f_layer)


def single_wedge(source):
    # note
    # it does not really make sense to do a two-hour interval
    # with a full sky,
    # but I do not want to risk large file sizes,
    # so the simulation will be deliberately truncated...
    
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    s_axis = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0]))
   
    for outer_ant in outer_ants.keys():
        inner_ants = outer_ants.copy()
        del inner_ants[outer_ant]

        for inner_ant in inner_ants.keys():
        
            f_layer = []
            for ni in nu_rl:
                s = s_axis[ni]
                A_n = A_full[ni]

                t_layer = []
                for ti in t_rl:
                    t = t_axis[ti]

                    r = rot.radec2lm(ra, dec, ra0=t)
                    phi = ant.phase_factor(
                        outer_ant, inner_ant, r, nu)
                    
                    A_t = A_n[ti]
                    
                    next_vista = np.dot(np.dot(A_t, s), phi)
                    t_layer.append(next_vista)
          
                f_layer.append(np.array(t_layer))

            inner_ants[inner_ant] = np.array(f_layer)

        outer_ants[outer_ant] = inner_ants

    return outer_ants
    
def package(block, ptitle):
    """
    Returns a print-ready dictionary,
        particularly for use with the pps.py routines.
        
    Use pickle_dict to save to disk.
    """
    return {'frequencies' : nu_axis,
            'times' : t_axis,
            'picture' : block,
            'title' : ptitle}

def pickle_dict(dict_, label):
    """
    Pickles @dict_ with a file name based on @label.
    While this is a generic routine and will pickle any dictionary,
        the intent is solely for use in conjunction with the
        package_wedge routine.
    """
    with open(label + '.pickle', 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
def auto_wedge(list_sources, label, ptitle):
    """
    Automatically runs
        full_wedge(list_sources)
        package_wedge()
        pickle_dict
    """
    pickle_dict(package(
        full_wedge(list_sources), ptitle
    ), label)
    
def auto_helix(ant1, ant2, sources, label, ptitle):
    pickle_dict(package(
        multi_helix(ant1, ant2, sources), ptitle
    ), label)
    
