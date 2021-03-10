import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

def auto_show(fname, static=False):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results using a
        2D + color wedge plot.
    """
    fcd, fs, ts, ptitle = load_wedge_sim(fname + ".pickle")
    print("Simulation file loaded.\n")
    
    transformed = transform_wedge(fcd, fs, ts)
    print("Fourier transforms applied to simulation.\n")
    
    if static:
        wedge = static_visual(sim_dict)
    else:
        wedge = collect_wedge_points(transformed, fs, ts)
    print("Wedge points collected.\n")
    
    plot_3D(wedge, ptitle)

def build_fourier_candidates(fname):
    sim_file = open(fname + ".pickle", "rb")
    
    meta = pickle.load(sim_file)
    
    ptitle = meta['title']
    
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
    
    return fourierc, raw_vis, fs, ts, ptitle
   
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
                np.log10(norm)
            )))
            
    return np.array(visual)    
    
def show_helix(fname, pt_override=None):
    """
    Load simulated visibilities from the file named
        @fname
    and visually interpret the results as a
    delay-spectrum helix a la (Parsons, 2012).
    
    @pt_override is a way to override the title with which
    a simulation came. It is extremely bad form to use this,
    but it can come in handy during lapses of diligence.
    """
    fourierc, raw_vis, fs, ts, ptitle = \
        build_fourier_candidates(fname)
    
    if pt_override is not None:
        ptitle = pt_override
    
    print("Data from file re-organized.")
    
    fouriered = transform_power(fourierc, fs, ts)
    
    print("Fourier transforms computed.")
    
    visual = collect_helix_points(fouriered, fs, ts)  
   
    print("Points collected.")
        
    plt.title("88m, 200 MHz bandwidth")
    plt.xlabel("Delays [ns]")
    plt.ylabel("LST [hr]")
    
    plot_3D(visual, ptitle)
    
    return visual, fouriered, fs, ts, raw_vis

def plot_3D(visual, title, scaled=False):
    """
    Primitive 3D plotter.
    
    For use with the return value of either
        static_wedge_vis
    or
        dynamic_wedge_vis
        
    Disable @scaled if you are using values such as logarithms
    """
    x = visual[:, 0]
    y = visual[:, 1]
    z = visual[:, 2]

    colors = None
    if (scaled):
        scaled_z = (z - z.min()) / z.ptp()
        colors = plt.cm.viridis(scaled_z)
    else:
        colors = z

    plt.title(title)

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

    ptitle = meta['title']

    fs = meta['frequencies']
    num_f = len(fs)
    
    ts = meta['times']
    num_t = len(ts)
    
    sim = meta['picture']
    
    fcd = {} # fourier candidate dictionary
    
    for ant1 in sim.keys():
        fcd[ant1] = {}
        for ant2 in sim[ant1].keys():
            # 0: I    1: Q    2: U    3: V
            fourierc = [[], [], [], []]
            
            for ti in range(num_t):
                for parameter in fourierc:
                    parameter.append([])
                
                for ni in range(num_f):
                    v = sim[ant1][ant2][ni][ti]

                    for p_idx in range(len(fourierc)):
                        fourierc[p_idx][ti].append(v[p_idx])
                
                for parameter in fourierc:
                    parameter[ti] = np.array(parameter[ti])
        
            for parameter in fourierc:
                parameter = np.array(parameter)

            fcd[ant1][ant2] = np.array(fourierc)

    return fcd, fs, ts, ptitle

def static_wedge_vis(sim_dict, fs):
    raise NotImplementedError("I have only updated" + \
    " the dynamic wedge routine at this moment.")
    
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
    
    etas = f2etas(fs)    

    center_f = np.average(fs)
    z = fq2z(center_f)
    lambda_ = C / center_f

    k_par = k_parallel(etas, z)
    k_starter = k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis
    
    ### aliasing ###
    nu_idxs = sim_dict['fs']

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
def collect_wedge_points(fcd, fs, ts):
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
    num_t = len(ts)
    num_f = len(fs)
    
    visual = []
    
    etas = pol.f2etas(fs)    

    center_f = np.average(fs)
    z = pol.fq2z(center_f)
    lambda_ = pol.C / center_f

    k_par = pol.k_parallel(etas, z)
    k_starter = pol.k_perp(z) / lambda_ # this will need to be
    # multiplied on a per-baseline basis

    for ant1 in fcd.keys():
        for ant2 in fcd[ant1].keys():
            k_orth = k_starter * sf.ant.baselength(ant1, ant2)
            
            for nu_idx in range(num_f):
                ###!!! Just using I at the moment
                powers_prop = []
                
                for t_idx in range(num_t - 1):
                    this_instant = \
                        fcd[ant1][ant2][:, t_idx, nu_idx]
                    next_instant = \
                        fcd[ant1][ant2][:, t_idx + 1, nu_idx]
                    
                    # this is just a proportionality.
                    powers_prop.append(np.abs(np.vdot(
                        this_instant,
                        next_instant
                    )))

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
                
                """
                dspecvec = np.array([
                parameter[ti][ni] for parameter in fouriered
                ])
            
                norm = np.linalg.norm(dspecvec)

                visual.append(np.array((
                    etas[ni] * 1e9,
                    ts[ti] * 12 / np.pi,
                    np.log10(norm)
                )))
                """

                visual.append(wedge_datum)   

    visual = np.array(visual)
   
    return np.array(visual)

### This is a really bad ad-hoc testing script.
### We want to scrap this ASAP
def micro_wedge(h1, f1, b1, h2, f2, b2, h3, f3, b3):
    """
    The axes do not line up with Nunhokee et al.
    Probably something wrong with your constants
        or usage thereof.
    """
    center_f1 = np.average(f1)
    z1 = pol.fq2z(center_f1)
    lambda1 = pol.C / center_f1
    k_par1 = pol.k_parallel(h1[:, 0], z1)
    k_orth1 = pol.k_perp(z1) / lambda1 * b1
    
    center_f2 = np.average(f2)
    z2 = pol.fq2z(center_f2)
    lambda2 = pol.C / center_f2
    k_par2 = pol.k_parallel(h2[:, 0], z2)
    k_orth2 = pol.k_perp(z2) / lambda2 * b2
    
    center_f3 = np.average(f3)
    z3 = pol.fq2z(center_f3)
    lambda3 = pol.C / center_f3
    k_par3 = pol.k_parallel(h3[:, 0], z3)
    k_orth3 = pol.k_perp(z3) / lambda3 * b3
    
    y = np.concatenate((k_par1, k_par2, k_par3))
    x = np.concatenate((
        np.repeat(k_orth1, len(k_par1)),
        np.repeat(k_orth2, len(k_par2)),
        np.repeat(k_orth3, len(k_par3))
    ))
    
    colors = np.concatenate((h1[:, 2], h2[:, 2], h3[:, 2]))

    plt.title("Helix concatenation")

    #print("Minimum:", z.min())
    #print("PTP:", z.ptp())

    plt.scatter(x, y, marker='.', c=colors)
    plt.colorbar()
    plt.show()
    
def transform_wedge(original, fs, ts):
    num_f = len(fs)
    num_t = len(ts)
    
    import copy
    
    fourier_dict = copy.deepcopy(original)

    window = pol.genWindow(num_f)

    for ant1 in fourier_dict.keys():
        for ant2 in fourier_dict[ant1].keys():
            fourier = fourier_dict[ant1][ant2]
            for ti in range(num_t):
                for parameter in fourier:
                    parameter[ti] = np.fft.fft(parameter[ti] * window)
                
    return fourier_dict
    
