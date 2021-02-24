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
    show_wedge(wedge)
    
def show_helix(fname):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results as a
    delay-spectrum helix a la (Parsons, 2012).
    """
    sim_file = open(fname + ".pickle", "rb")
    
    meta = pickle.load(sim_file)
    
    frq = meta['frequencies']
    etas = f2etas(frq)
    
    ts = meta['times']
    #ts = np.fft.fft(meta['times'])
    #plt.plot(ts); plt.show()
    
    sim = meta['picture']
    
    # 0: I    1: Q    2: U    3: V
    fourier = [[], [], [], []]
    
    raw_vis = []
    
    for ti in range(len(sim[0])):
        for parameter in fourier:
            parameter.append([])
        
        raw_vis.append([])
        
        last_i = len(fourier[0]) - 1
        
        for ni in range(len(sim)):
            
            v = sim[ni][ti]
            # print(v)
            
            for p_idx in range(len(fourier)):
                fourier[p_idx][last_i].append(v[p_idx])
            
            raw_vis[last_i].append(sim[ni][ti])
            #norm = np.linalg.norm(sim[ni][ti])
        
        for parameter in fourier:
            parameter[last_i] = np.array(parameter[last_i])
        
        raw_vis[last_i] = np.array(raw_vis[last_i])
        
    for parameter in fourier:
        parameter = np.array(parameter)

    fourier = np.array(fourier)
    
    raw_vis = np.array(raw_vis)
    
    print("Data from file re-organized.")
    
    N = len(frq)
    df = frq[1] - frq[0]
    window = genWindow(N)
        
    print("Window generated.")
    
    visual = []
    
    # option 10
    # ts = np.fft.fftshift(ts)
    
    for ti in range(int(len(fourier[0]) / 2)):
        #print(fourier_i[ti])
        #print("Length is", len(fourier_i[ti]))
        
        """
        # option 6 [next p lines]
        fourier_i[ti] = \
            np.fft.fftshift(np.fft.fft(fourier_i[ti]))
        fourier_q[ti] = \
            np.fft.fftshift(np.fft.fft(fourier_q[ti]))
        fourier_u[ti] = \
            np.fft.fftshift(np.fft.fft(fourier_u[ti]))
        fourier_v[ti] = \
            np.fft.fftshift(np.fft.fft(fourier_v[ti]))
        """
        
        """
        # what I had been doing before 2/17/21
        # aka option 5 [next 4 lines]
        fourier_i[ti] = np.fft.fft(fourier_i[ti])
        fourier_q[ti] = np.fft.fft(fourier_q[ti])
        fourier_u[ti] = np.fft.fft(fourier_u[ti])
        fourier_v[ti] = np.fft.fft(fourier_v[ti])
        """
        
        # fft with window: option 9 [next 4 lines]
        for parameter in fourier:
            parameter[ti] = np.fft.fft(parameter[ti] * window)
        
        """
        # ifft: option 7 [next p lines]
        fourier_i = np.fft.ifft(fourier_i)
        fourier_q = np.fft.ifft(fourier_q)
        fourier_u = np.fft.ifft(fourier_u)
        fourier_v = np.fft.ifft(fourier_v)
        """
        
        """
        # ifft with window: option 8 [next 4 lines]
        fourier_i = np.fft.ifft(fourier_i * window)
        fourier_q = np.fft.ifft(fourier_q * window)
        fourier_u = np.fft.ifft(fourier_u * window)
        fourier_v = np.fft.ifft(fourier_v * window)
        """
        
        percent = (ti + 1) * 100 / len(fourier[0])
        print("Fourier transforms", percent, "% complete.")
        
        #print("Fourier transforms computed")
        
        for ni in range(len(fourier[0][ti])):
            dspecvec = np.array([
                parameter[ti][ni] for parameter in fourier
            ])
        
            alpha_norm = np.linalg.norm(dspecvec)
            beta_norm = float(np.dot(np.conj(dspecvec), dspecvec))
            #print(beta_norm)
        
            #!! Are all the indices lining up as I think they are?
            visual.append(np.array((
                etas[ni], ts[ti],
                alpha_norm
            )))
        
    visual = np.array(visual)
    # visual = np.fft.fftshift(visual)
    
    delays = visual[:, 0] * 1e9
    #delays = np.fft.fftshift(delays)
    
    times = visual[:, 1] * 12 / np.pi
    #times = np.fft.fftshift(times)
    
    v = visual[:, 2]
    v = np.log10(v)
    #v = np.fft.fftshift(v)
    
    print("t zero", times[0])

    #scaled_v = (v - v.min()) / v.ptp()
    #colors = plt.cm.viridis(scaled_v)
    colors = plt.cm.viridis(v)

    plt.title("88m, 200 MHz bandwidth")

    " We HAVE to do better than this. How do I line up a color bar? "

    print("Minimum:", v.min())
    print("PTP:", v.ptp())

    plt.xlabel("Delays [ns]")
    plt.ylabel("LST [hr]")

    #print(len(delays))
    
    #plt.plot(delays[:len(frq)], v[:len(frq)])
    #plt.show()
    
    plt.scatter(delays, times, marker='.', c=colors)
    
    plt.colorbar()
    plt.show()
    
    return fourier, raw_vis
    #return delays[:len(frq)], raw_vis[:len(frq)], fourier_field
    
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
    
def show_wedge(wedge):
    """
    Primitive 3D plotter.
    
    For use with the return value of either
        static_wedge_vis
    or
        dynamic_wedge_vis
    """
    k_orth = wedge[:, 0]
    k_parr = wedge[:, 1]
    p_p = wedge[:, 2]

    scaled_pow = (p_p - p_p.min()) / p_p.ptp()
    colors = plt.cm.viridis(scaled_pow)

    print("Minimum:", p_p.min())
    print("PTP:", p_p.ptp())

    plt.scatter(k_orth, k_parr, marker='.', c=colors)
    plt.colorbar()
    plt.show()

