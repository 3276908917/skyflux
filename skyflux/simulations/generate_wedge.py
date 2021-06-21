import matplotlib.pyplot as plt
import numpy as np
import pickle

from skyflux import catalog
from skyflux import ant
from skyflux import vis
from skyflux import stokes
from skyflux import rot
from skyflux import demo
from skyflux import utils

# disgustingly hacky
MACRO_EPSILON = 0.001

# constants
HOUR = 2 * np.pi / 24
MINUTE = HOUR / 60
SECOND = MINUTE / 60

# we keep these as global parameters to avoid the potential overhead
# of passing by value

# For wedges
# 8 MHz from 145 to 155, split them into one MHz resolution
nu_axis = np.arange(125e6, 175e6 + MACRO_EPSILON, 1e6)
t_axis = np.arange(3 * HOUR, 5 * HOUR, 4 * MINUTE)
    # it's an arbitrary region of the cold patch

nu_rl = range(len(nu_axis))
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
        J_source = stokes.create_J(
            az=azs, alt=alts, nu=nu, radians=True)
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
    Contains every possible frequency A matrix;
    When performing calculations at frequency index f, we say
        this_A = A_tensor[f]
    """
    A_tensor = []

    az, alt = rot.eq_to_topo(ra, dec, lst=lst0, radians=True)

    for nu in nu_axis:
        J_source = stokes.create_J(
            az=az, alt=alt, nu=nu, radians=True)
        A_source = stokes.create_A(J=J_source)

        A_tensor.append(np.array(A_source))
        
    return np.array(A_tensor)

# Scan over all frequencies, for a single source, over all possible baselines
#! Is this method still relevant?
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
                   
def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nSimulation: " + percent_status + "% complete.")

def null_source(obj):
    return obj.alpha != obj.alpha

def full_wedge(sources=catalog.srcs):
    percent_interval = 100 / len(sources)
    percent = 0
    
    wedge = single_wedge(sources[0])
    
    percent += percent_interval
    tick(percent)
    
    for next_obj in sources[1:]:
        if null_source(next_obj):
            continue
        
        next_wedge = single_wedge(next_obj)
        wedge = merge_wedges(wedge, next_wedge)

        percent += percent_interval
        tick(percent)
        
    return wedge

outer_ants = ant.ant_pos.copy()

def single_wedge(source):
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
                nu = nu_axis[ni]
            
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
          
                t_layer = np.array(t_layer)
                
                f_layer.append(t_layer)

            inner_ants[inner_ant] = np.array(f_layer)

            f_layer = np.array(f_layer)
            
        outer_ants[outer_ant] = inner_ants

    return outer_ants
    
def single_wedge_readout(source):
    """
    Return the time indices at which
    the source is most visible
    over the first three hypothetical plots.
    
    Chiefly for use with the function
    single_wedge_spot 
    """
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    s_axis = []

    kill_timer = 0
    
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
                nu = nu_axis[ni]
            
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
          
                t_layer = np.array(t_layer)
                
                f_layer.append(t_layer)

            inner_ants[inner_ant] = np.array(f_layer)

            f_layer = np.array(f_layer)
            
            max_ti = 0
            max_V = float('-inf')
            
            for ti in range(len(f_layer[0])):
                if f_layer[0, ti, 0] > max_V:
                    max_ti = ti
                    max_V = f_layer[0, ti, 0]
        
            print(max_ti)                    
            
            kill_timer += 1
            
            if kill_timer > 3:
                return
            
        outer_ants[outer_ant] = inner_ants

    return outer_ants

def single_wedge_spot(source, ti):
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_full = A_tensor(ra, dec)

    s_axis = []
    
    for ni in nu_rl:
        nu = nu_axis[ni]
        I = vis.get_I(source, nu)
        s_axis.append(np.array([complex(I), 0, 0, 0]))

    kill_timer = 0
   
    for outer_ant in outer_ants.keys():
        inner_ants = outer_ants.copy()
        del inner_ants[outer_ant]

        for inner_ant in inner_ants.keys():
        
            f_layer = []
            for ni in nu_rl:
                nu = nu_axis[ni]
            
                s = s_axis[ni]
                A_n = A_full[ni]

                t_layer = []
                
                t = t_axis[ti[kill_timer]]

                r = rot.radec2lm(ra, dec, ra0=t)
                phi = ant.phase_factor(
                    outer_ant, inner_ant, r, nu)
                
                A_t = A_n[ti[kill_timer]]
                
                Jonesed = np.dot(A_t, s)
                
                next_vista = np.dot(Jonesed, phi)
                
                f_layer.append(np.array([Jonesed, next_vista]))

            f_layer = np.array(f_layer)

            inner_ants[inner_ant] = f_layer
                
            """
            * Check the s input vector
            * Multiply s by the Jones matrix
            * Output visibilities
            """
            
            ad_hoc_labels = ["I", "Q", "U", "V"]
            
            log_enabled = False
            
            x = nu_axis / 1e6
            
            s_axis = np.array(s_axis)
            
            if log_enabled:
                s_plot = np.log10(s_axis[:, 0])
            else:
                s_plot = s_axis[:, 0]
            
            # Jonesed
            for i in range(4):
                y = f_layer[:, 0, i]
                
                if log_enabled:
                    plt.plot(
                        x,
                        np.log10(np.abs(y)),
                        label=ad_hoc_labels[i]
                    )
                else:
                    plt.plot(
                        x,
                        np.abs(y),
                        label=ad_hoc_labels[i]
                    )
            
            """
            plt.plot(
                x,
                s_plot,
                label="s vector, I"
            )
            """
            
            plt.legend(loc='upper right')
            plt.title(
                "Wedge generator, no phase factor. Antennae: " + \
                str(outer_ant) + " to " + \
                str(inner_ant)
            )
            plt.xlabel("Frequency [MHz]")
            plt.ylabel("Brightness Magnitude [Jy]")
            plt.show()
            
            # Jonesed
            for i in range(4):
                y = f_layer[:, 1, i]
                
                if log_enabled:
                    plt.plot(
                        x,
                        np.log10(np.abs(y)),
                        label=ad_hoc_labels[i]
                    )
                else:
                    plt.plot(
                        x,
                        np.abs(y),
                        label=ad_hoc_labels[i]
                    )
            
            """
            * plot A at the point of interest
                * do Q/I, U/I ratios 
            * plot A multiplied by s
            
            """
            
            """
            plt.plot(
                x,
                s_plot,
                label="s vector, I"
            )
            """
            
            plt.legend(loc='upper right')
            plt.title(
                "Wedge generator, with phase factor. Antennae: " + \
                str(outer_ant) + " to " + \
                str(inner_ant)
            )
            plt.xlabel("Frequency [MHz]")
            plt.ylabel("Brightness Magnitude [Jy]")
            plt.show()
            
            kill_timer += 1
            
            if kill_timer > 3:
                return

        outer_ants[outer_ant] = inner_ants

    return outer_ants
    
def merge_wedges(wedge1, wedge2):
    """ We assume that both wedges have the same format:
        ant1, ant2, nu, t hierarchies are exactly the same.
    """
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
                    
                    sum_[ant1][ant2][nu_idx][t_idx] = \
                        np.add(vis1, vis2)
                        
    return sum_
    
def merge_files(fname1, fname2, new_fname, new_ptitle):
    """ This function only really makes sense if f1 and f2 
        have common simulation parameters
            (frequency res, time res, etc)
    """
    f1 = open(fname1 + ".pickle", "rb")
    f2 = open(fname2 + ".pickle", "rb")
    
    sim1 = pickle.load(f1)
    sim2 = pickle.load(f2)
    
    # arbitrarily take common sim. parameters from f1
    fs = sim1['frequencies']
    ts = sim1['times']
    
    # merge the wedges
    wedge1 = sim1['picture']
    wedge2 = sim2['picture']
    wmerged = merge_wedges(wedge1, wedge2)
    
    packed = {
        'frequencies' : fs,
        'times' : ts,
        'picture' : wmerged,
        'title' : new_ptitle
    }
    pickle_dict(packed, new_fname)
    
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
        
def auto_wedge(list_sources, label, ptitle):
    """
    Automatically runs
        full_wedge(list_sources)
        package_wedge()
        pickle_dict
    """
    utils.pickle_dict(package(
        full_wedge(list_sources), ptitle
    ), label)
    
