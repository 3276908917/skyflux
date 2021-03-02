#! /usr/bin/env python

""" Begin section: pared-down copy of
    Chuneeta/PolarizedSims/COSMO_constants.py """

import numpy as np
import astropy.cosmology as CS

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1.42040575177 # FREQUENCY OF 21 CM HYDROGEN LINE IN GHZ
COSMO = CS.FlatLambdaCDM(H0=100.0, Om0=0.27) # H0 = 100 km/s/Mpc

def fq2z(fq):
   """
   Redshift corresponding the specified frequency

   Input(s)
      fq :  [scalar] frequency in Hz
   """
   return F21 / fq - 1
   
def transverse_comoving_distance(z):
   """
   Transverse comoving distance at redshift z
       corresponding to an
       angular separation of 1 radian in Mpc/h.
   
   Input(s)
      z :  [scalar] redshift
   """
   return COSMO.comoving_distance(z).value # Mpc/h 

def comoving_depth(B,z):
   """
   Comoving line-of-sight depth corresponding to
       specified redshift and bandwidth for redshifted
       21 cm line in Mpc/h
  
   Input(s)
      B :    [scalar] Observing bandwith in Hz
      z :    [scalar] redshift
   """
   return (C/1e3) * B * (1+z)**2 / \
       F21/COSMO.H0.value/COSMO.efunc(z) # Mpc/h 

def dkprll_deta(z):
   """
   Constant to transform delays to line-of-sight
       wavenumbers corresponding to redshift and
       21 CM HI line in h/Mpc
   
   Input(s)
      z :  [scalar] redshift
   """
   return 2 * np.pi * COSMO.H0.value * F21 * \
       COSMO.efunc(z) / C /  (1+z)**2 * 1e3

def k_parallel(delays, z):
   """
   Compute line-of-sight wavenumbers corresponding to
       specified delays and redshift for
       redshifted 21 cm line in h/Mpc

   Input(s):
      z : [scalar] redshift
   """
   return dkprll_deta(z) * delays

def k_perp(z):
   """
   Compute transverse wavenumbers corresponding to
       redshifted 21 cm line in h/Mpc

   Input(s)
      z              : [scalar] redshift
   """
   return 2 * np.pi / transverse_comoving_distance(z)

""" End section: Chuneeta/PolarizedSims/COSMO_constants.py """

""" Begin section: pared-down copy of
    Chuneeta/PolarizedSims/genPowerSpectra.py"""

def genWindow(size):
    """
    Implements Blackmann-Harris filter

    size : Size/Lenth of frequency channel to which
        filter is applied; type:int
    """
    window = np.zeros((size),)
    alpha = 0.16 #! magic number
    _w = lambda n: \
        (1-alpha) / 2. - 0.5 * np.cos((2 * np.pi * n)/(size - 1)) \
        + alpha / 2. * np.cos((4 * np.pi * n)/(size - 1))
    window[:] = _w(np.arange(size))
    return window

def f2etas(freqs):
    """
    Evaluates geometric delay (fourier conjugate of frequency)
   
    -freqs: Frequencies in GHz; type:numpy.ndarray
    """
    df = freqs[1] - freqs[0]
    etas = np.fft.fftfreq(freqs.size, df)
    return etas

def delay_transform(data,fqs,convert=None):
    """
    Fourier transforms visibility along frequency axis

    - data: per baseline visibility; type:numpy.ndarray
    - fqs:  slected frequencies in GHz; dtypw:numpy.ndarray
    """
    N = fqs.size
    df = fqs[1] - fqs[0]
    window = genWindow(N) #! this label is never used
    delaySpec = np.fft.ifft(data) * N * df
    return delaySpec 
      
""" End section: Chuneeta/PolarizedSims/genPowerSpectra.py """

""" Remainder of this file: the code unique to this script.
    i.e. my own work. "power plot sketch" """

import skyflux as sf

import matplotlib.pyplot as plt
import pickle

def auto_show(fname, static=False):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results using a
        2D + color wedge plot.
    """
    sim_dict = load_wedge_sim(fname + ".pickle")
    if static:
        wedge = static_visual(sim_dict)
    else:
        wedge = dynamic_visual(sim_dict) 
    plot_3D(wedge)

def build_fourier_candidates(fname):
    sim_file = open(fname + ".pickle", "rb")
    
    meta = pickle.load(sim_file)
    
    fs = meta['frequencies']
    num_f = len(fs)
    
    ts = meta['times']
    num_t = len(ts)
    
    sim = meta['picture']
    
    # 0: I    1: Q    2: U    3: V
    fourierc = [[], [], [], []]
    
    raw_vis = []
    
    for ti in range(num_t):
        for parameter in fourierc:
            parameter.append([])
        
        raw_vis.append([])
        
        for ni in range(num_f):
            v = sim[ni][ti]

            for p_idx in range(len(fourierc)):
                fourierc[p_idx][ti].append(v[p_idx])
            
            raw_vis[ti].append(v)
            #norm = np.linalg.norm(sim[ni][ti]) same outcome
        
        for parameter in fourierc:
            parameter[ti] = np.array(parameter[ti])
        
        raw_vis[ti] = np.array(raw_vis[ti])
        
    for parameter in fourierc:
        parameter = np.array(parameter)

    fourierc = np.array(fourierc)
    
    raw_vis = np.array(raw_vis)
    
    return fourierc, raw_vis, fs, ts
   
def transform_power(original, fs, ts):
    num_f = len(fs)
    num_t = len(ts)
    
    import copy
    
    fourier = copy.deepcopy(original)

    window = genWindow(num_f)

    for ti in range(num_t):
        """
        # option 6
        for parameter in fourier:
            parameter[ti] = \
                np.fft.fftshift(np.fft.fft(parameter[ti])
        """
        
        """
        # what I had been doing before 2/17/21
        # aka option 5
        for parameter in fourier:
            parameter[ti] = np.fft.fft(parameter[ti])
        """
        
        # fft with window: option 9
        for parameter in fourier:
            parameter[ti] = np.fft.fft(parameter[ti] * window)
        
        """
        # ifft: option 7
        for parameter in fourier:
            parameter[ti] = np.fft.ifft(parameter[ti])
        """
        
        """
        # ifft with window: option 8 [next 4 lines]
        for parameter in fourier:
            parameter[ti] = np.fft.ifft(parameter[ti] * window)
        """
        
    return fourier 

def collect_helix_points(fouriered, fs, ts):
    num_t = len(ts)
    num_f = len(fs)
    
    visual = []
    
    etas = f2etas(fs)
        
    for ti in range(num_t):
        for ni in range(num_f):
            dspecvec = np.array([
                parameter[ti][ni] for parameter in fouriered
            ])
        
            norm = np.linalg.norm(dspecvec)

            visual.append(np.array((
                etas[ni] * 1e9,
                ts[ti] * 12 / np.pi,
                norm #np.log10(norm)
            )))
            
    return np.array(visual)    
    
def show_helix(fname):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results as a
    delay-spectrum helix a la (Parsons, 2012).
    """
    fourierc, raw_vis, fs, ts = build_fourier_candidates(fname)
    
    print("Data from file re-organized.")
    
    fouriered = transform_power(fourierc, fs, ts)
    
    print("Fourier transforms computed.")
    
    visual = collect_helix_points(fouriered, fs, ts)  
   
    print("Points collected.")
        
    plt.title("88m, 200 MHz bandwidth")
    plt.xlabel("Delays [ns]")
    plt.ylabel("LST [hr]")
    
    plot_3D(visual)
    
    return fouriered, fs, ts, raw_vis

def plot_3D(visual):
    """
    Primitive 3D plotter.
    
    For use with the return value of either
        static_wedge_vis
    or
        dynamic_wedge_vis
    """
    x = visual[:, 0]
    y = visual[:, 1]
    z = visual[:, 2]

    scaled_z = (z - z.min()) / z.ptp()
    colors = plt.cm.viridis(scaled_z)

    print("Minimum:", z.min())
    print("PTP:", z.ptp())

    plt.scatter(x, y, marker='.', c=colors)
    plt.colorbar()
    plt.show()
    
def load_wedge_sim(fname):
    """
    Fill out a wedge data structure based on
        the contents of the file with name @fname.
    The returned dictionary has the following key-value pairs:
        "fs": the frequencies used in the simulation
        "kp": array of k-parallel, the proxy-variable for frequency
        "ks": a constant coefficient, to be multiplied by
            baseline lengths to yield k-perpendicular,
            the proxy-variable for baseline
        "sim": a three-dimensional array of visibilities
            the first index corresponds to frequency
            the second index corresponds to LST
            the third index corresponds to Stokes parameter
                0: I
                1: Q
                2: U
                3: V
    """
    sim_file = open(fname, "rb")

    meta = pickle.load(sim_file)

    frq = meta['frequencies'] #! I don't remember where this
       # naming scheme is in the code

    # f2etas claims to accept frequencies in GHz,
    # but the axes come out strangely
    etas = f2etas(frq)# / 1e9)
    
    # we are ranging from 50 to 250 MHz
    center_f = np.average(frq)
    z = fq2z(center_f)
    lambda_ = C / center_f

    k_par = k_parallel(etas, z)
    k_starter = k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis

    nu_idxs = range(len(frq))

    return {"fs": nu_idxs,
            "kp": k_par,
            "ks": k_starter,
            "sim": meta['picture']}

def static_wedge_vis(sim_dict):
    """
    Read from the wedge data structure @sim_dict
        (specifically, one using the format
        established in load_wedge_sim)
    and generate 3D points appropriate for a wedge plot.
    i.e., return a list of triples:
        (k_perpendicular, k_parallel, power*)
    
    * currently, the implementation uses values that
    should be proportional to power. The final constants
    have not yet been considered.
    
    This function is distinct from
        dynamic_wedge_vis
    in assuming that the simulation corresponding to
        @sim_dict
    runs over only one value of LST.
    """
    wedge_data = []
    
    ### aliasing ###
    nu_idxs = sim_dict['fs']
    k_par = sim_dict['kp']
    k_starter = sim_dict['ks']
    sim = sim_dict['sim']

    for ant1 in sim.keys():
        for ant2 in sim[ant1].keys():
            """
            Since we only used one LST,
            we do not need to do any averaging
                (averaging is supposed to happen over LSTs,
                    not frequency)
            """
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            for nu_idx in nu_idxs:
                # this is a proportionality.
                # The real deal uses the power equation 6
                    # from Nunhokee et al.
                    
                brightness = sim[ant1][ant2][nu_idx]
                    
                power_prop = np.log10(np.vdot(
                    brightness,
                    brightness
                ))
                
                wedge_datum = np.array([
                    k_orth,
                    k_par[nu_idx],
                    float(power_prop)
                ])

                wedge_data.append(wedge_datum)
                
    return np.array(wedge_data)
    
#! There has to be some way to merge this with the function
# above
def dynamic_wedge_vis(sim_dict):
    """
    Read from the wedge data structure @sim_dict
        (specifically, one using the format
        established in load_wedge_sim)
    and generate 3D points appropriate for a wedge plot.
    i.e., return a list of triples:
        (k_perpendicular, k_parallel, power*)
    
    * currently, the implementation uses values that
    should be proportional to power. The final constants
    have not yet been considered.
    
    This function is distinct from
        static_wedge_vis
    in accepting multiple LST values from
        @sim_dict.
    """
    wedge_data = []
    
    ### aliasing ###
    nu_idxs = sim_dict['fs']
    k_par = sim_dict['kp']
    k_starter = sim_dict['ks']
    sim = sim_dict['sim']

    for ant1 in sim.keys():
        for ant2 in sim[ant1].keys():
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            
            for nu_idx in nu_idxs:
                system = sim[ant1][ant2][nu_idx]
                powers_prop = []
                
                for t_idx in range(len(system) - 1):
                    this_instant = system[t_idx]
                    next_instant = system[t_idx + 1]
                    # this is a proportionality.
                    # The real deal uses the power equation 6
                    # from Nunhokee et al.
                    powers_prop.append(np.vdot(
                        this_instant,
                        next_instant
                    ))

                    # Does v have units of brightness or square of brightness?
                    # We want a delay transform delay_transform(data, fqs)

                wedge_datum = np.array([
                    k_orth,
                    k_par[nu_idx],
                    float(
                        np.log10(
                            np.average(
                                np.array(powers_prop)
                            )
                        )
                    )
                ])

                wedge_data.append(wedge_datum)   

    ### begin prototyping section

    wedge_data = np.array(wedge_data)

    """
       The structure of wedge_data is a collection of triplets
       (k_orth, k_parr, some squared brightness)

       Based on the order in which I added data, I think that each consecutive
       len(nu_idxs) should have the same value of k_orth, right?
    """
    
    """
    nu_rl = len(nu_idxs)

    wedge_final = wedge_data.copy()

    for i in range(len(wedge_data) / nu_rl - 1):
        start = i * nu_rl
        end = (i + 1) * nu_rl
        vis_over_f = wedge_data[start:end]
        dt = delay_transform(vis_over_f[:, 2], frq / 1e9)
        # now copy that into wedge_final somehow
        for j in range(len(vis_over_f)):
            wedge_final[j + start][2] = vis_over_f[j]
       
    ### terminate prototyping section  
    """
   
    return np.array(wedge_data)
    
