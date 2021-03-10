import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

def hauto_show(fname, pt_override=None):
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

    window = pol.genWindow(num_f)

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
    
    etas = pol.f2etas(fs)
        
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
    
