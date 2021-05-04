import numpy as np
import matplotlib.pyplot as plt
import pickle

import skyflux as sf
import skyflux.deprecated.polSims as pol

def probe_stokes_params(fname, sp=None):
    """
    out of context
    """
    fcd, fs, ts, ptitle = load_wedge_sim(fname + ".pickle")
    
    print("Simulation file loaded.\n")
    
    transformed = transform_wedge(fcd, fs, ts)
    print("Fourier transforms applied to simulation.\n")
    
    wedge = collect_wedge_points(
        transformed, fs, ts, sp
    )
    
def load_wedge_sim(fname):
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
            fourierc = [[], [], [], []]
            
            for ti in range(num_t):
                for parameter in fourierc:
                    parameter.append([])
                
                for ni in range(num_f):
                    v = sim[ant1][ant2][ni][ti]

                    v = np.array([
                        np.abs(S) for S in v
                    ])
                    
                    for i in range(len(v)):
                        if v[i] == v.max():
                           max_counts[i] += 1 

                    for p_idx in range(len(fourierc)):
                        fourierc[p_idx][ti].append(v[p_idx])
                
                for parameter in fourierc:
                    parameter[ti] = np.array(parameter[ti])
        
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
                print(len(fourier))
                for parameter in fourier:
                    parameter[ti] = np.fft.fft(
                        parameter[ti] * window
                    )
                
    return fourier_dict

def collect_wedge_points(fcd, fs, ts, sp=None):
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

    for ant1 in fcd.keys():
        for ant2 in fcd[ant1].keys():
        
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
                    
                    this_instant = np.array([
                        np.abs(S) for S in this_instant
                    ])
                    
                    for i in range(len(this_instant)):
                        if this_instant[i] == this_instant.max():
                           max_counts[i] += 1 
                    
    print(max_counts)
