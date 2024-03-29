import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

ad_hoc_labels = ["I", "Q", "U", "V"]
ad_hoc_colors = ['k', 'b', 'orange', 'g']

def probe_stokes_params(fname, rAnt1=None, rAnt2=None):
    """
    out of context
    """
    fcd, fs, ts, ptitle = load_wedge_sim(
        fname + ".pickle", rAnt1, rAnt2
    )
    
    print("Simulation file loaded.\n")
    
    transformed = transform_wedge(fcd, fs, ts)
    print("Fourier transforms applied to simulation.\n")
    
    wedge = collect_wedge_points(
        transformed, fs, ts, rAnt1, rAnt2
    )
    
def load_wedge_sim(fname, rAnt1, rAnt2, tGoal=None):
    """
    out of context
    """
    sim_file = open(fname, "rb")

    meta = pickle.load(sim_file)

    ptitle = meta['title']

    fs = meta['frequencies']
    #print(fs)
    
    num_f = len(fs)
    
    ts = meta['times']
    num_t = len(ts)
    
    sim = meta['picture']
    
    fcd = {} # fourier candidate dictionary
    
    max_counts = [0, 0, 0, 0]
    
    for ant1 in sim.keys():
        fcd[ant1] = {}
        for ant2 in sim[ant1].keys():
            # 0: I    1: Q    2: U    3: V
            if rAnt1 is not None and rAnt1 != ant1:
                continue
            if rAnt2 is not None and rAnt2 != ant2:
                continue
            
            fourierc = [[], [], [], []]
            
            for ti in range(num_t):
                if tGoal is None or tGoal != np.around(ts[ti], 4):
                    continue
                
                for parameter in fourierc:
                    parameter.append([])
                
                for ni in range(num_f):
                    v = sim[ant1][ant2][ni][ti]

                    v_test = np.array([
                        np.abs(S) for S in v
                    ])
                    
                    for i in range(len(v_test)):
                        if v_test[i] == v_test.max():
                           max_counts[i] += 1 

                    for p_idx in range(len(fourierc)):
                        fourierc[p_idx][ti].append(v[p_idx])
                
                # si: Stokes index
                for si in range(len(fourierc)):
                    parameter = fourierc[si]
                    parameter[ti] = np.array(parameter[ti])
                    
                    # cut the next three statements
                    # in case you no longer want to plot this
                    plt.plot(fs / 1e6,
                        np.abs(parameter[len(parameter) - 1]),
                        label=str(ad_hoc_labels[si]),
                        color=ad_hoc_colors[si])
        
                plt.legend(loc='upper right')
                plt.title(
                    "Raw Visibilities (LST = " + \
                    str(np.around(ts[ti], 4)) + \
                    "). Antennae: " + \
                    str(ant1) + " to " + \
                    str(ant2)
                )
                plt.xlabel("Frequency [MHz]")
                plt.ylabel("Brightness Magnitude [Jy]")
                plt.show()
        
            for parameter in fourierc:
                parameter = np.array(parameter)

            fcd[ant1][ant2] = np.array(fourierc)

    print(max_counts)

    return fcd, fs, ts, ptitle    
       
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
                    parameter[ti] = np.fft.fft(
                        parameter[ti] * window
                    )
                
    return fourier_dict

def collect_wedge_points(fcd, fs, ts, rAnt1, rAnt2):
    """
    out of context
    """
    num_t = len(ts)
    num_f = len(fs)
    B = fs.max() - fs.min()
    
    visual = []
    
    etas = pol.f2etas(fs)
   
    k_para = np.array([
        pol.k_parallel(etas[i], pol.fq2z(fs[i] / 1e9)) \
        for i in range(num_f)
    ])

    max_counts = [0, 0, 0, 0]
    hmax_counts = [0, 0, 0, 0]

    for ant1 in fcd.keys():
        for ant2 in fcd[ant1].keys():
        
            if rAnt1 is not None and rAnt1 != ant1:
                continue
            if rAnt2 is not None and rAnt2 != ant2:
                continue
        
            baselength = sf.ant.baselength(ant1, ant2)
            
            for nu_idx in range(num_f):
                
                nua = np.average(fs)
                nu = fs[nu_idx]
                
                z = pol.fq2z(nu / 1e9)
                
                lambda_ = pol.C / nu
                
                k_perp = baselength * pol.k_perp(z) / lambda_
                
                for t_idx in range(num_t - 1):
                    this_instant = \
                        fcd[ant1][ant2][:, t_idx, nu_idx]
                    
                    # [I1, Q1, U1, V1] * [I2*, Q2*, U2*, V2*]
                    
                    this_test = np.array([
                        np.abs(S) for S in this_instant
                    ])
                    
                    for i in range(len(this_test)):
                        if this_test[i] == this_test.max():
                           max_counts[i] += 1 
                           
                    next_instant = \
                        fcd[ant1][ant2][:, t_idx + 1, nu_idx]
                    
                    # [I1, Q1, U1, V1] * [I2*, Q2*, U2*, V2*]
                    
                    hadamard = np.multiply(
                        this_instant, next_instant
                    )
                    
                    htest = np.array([
                        np.abs(S) for S in hadamard
                    ])
                    
                    for i in range(len(htest)):
                        if htest[i] == htest.max():
                           hmax_counts[i] += 1
                    
    print(max_counts)
    print(hmax_counts)
    
"""
Create a cross-sectional plot
    1. Just the visibilities, with frequency as axis
    2. Hadamard product of visibilities
    3. Multiply by the cosmological constant
    4. Divide by the normalization
"""

