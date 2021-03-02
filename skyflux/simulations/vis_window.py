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
# t_axis = np.arange(lst0 - hour, lst0 + hour, 4 * minute)
# For helices: 30 second intervals
t_axis = np.arange(0, 24 * HOUR, 30 * SECOND)
t_rl = range(len(t_axis))

def A_tensor(ra, dec, nu):
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

    J_source = stokes.create_J(az=azs, alt=alts, nu=nu, radians=True)
    A_source = np.array([stokes.create_A(J=J) for J in J_source])
        
    return np.array(A_source)
                   
def tick(percent):
    """ Give the user a progress update."""
    percent_status = str(np.around(percent, 4))
    print("\nSimulation: " + percent_status + "% complete.")

def null_source(obj):
    return obj.alpha != obj.alpha

def fwhm(ant1, ant2, source, nu, max_):
    """
    Full width at half maximum 
    """
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_n = A_tensor(ra, dec, nu)

    s_axis = []
        
    f_layer = []
        
    I = vis.get_I(source, nu)
    s = np.array([complex(I), 0, 0, 0])

    t_layer = []
    
    vnorm1 = None
    vnorm2 = None
    
    for ti in t_rl:
        t = t_axis[ti]

        A_t = A_n[ti]
        
        r = rot.radec2lm(ra, dec, ra0=t)
        phi = ant.phase_factor(ant1, ant2, r, nu)
        
        next_vista = np.dot(np.dot(A_t, s), phi)
        
        vis_norm = np.linalg.norm(next_vista)
        
        temp = vnorm2
        vnorm2 = vis_norm
        vnorm1 = temp
        
        if vnorm1 is None:
            continue
        elif 2 * vnorm2 >= max_ and 2 * vnorm1 <= max_:
            print("Rise time:", t, "radians")
            print("(" + rad_to_time(t) + ")\n")
        elif 2 * vnorm2 <= max_ and 2 * vnorm1 >= max_:
            print("Set time:", t)
            print("(" + rad_to_time(t) + ")")

def rad_to_time(theta):
    hr = int(theta / HOUR)
    rem = theta - HOUR * hr
    m = int(rem / MINUTE)
    rem = rem - MINUTE * m
    s = int(rem / SECOND)
    return str(hr) + " hours, " + \
        str(m) + " minutes, " + \
        str(s) + " seconds."

def vmax(ant1, ant2, source, nu):
    """
    Maximum visibility. For use with fwhm
    """
    max_ = float("-inf")
    
    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A_n = A_tensor(ra, dec, nu)

    s_axis = []
        
    f_layer = []
        
    I = vis.get_I(source, nu)
    s = np.array([complex(I), 0, 0, 0])

    t_layer = []
    
    for ti in t_rl:
        t = t_axis[ti]

        A_t = A_n[ti]
        
        r = rot.radec2lm(ra, dec, ra0=t)
        phi = ant.phase_factor(ant1, ant2, r, nu)
        
        next_vista = np.dot(np.dot(A_t, s), phi)
        
        vis_norm = np.linalg.norm(next_vista)
        
        if (vis_norm > max_):
            max_ = vis_norm
    
    return max_
    
def find_window(ant1, ant2, source, nu):
    m = vmax(ant1, ant2, source, nu)
    fwhm(ant1, ant2, source, nu, m)
    
