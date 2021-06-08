"""
TODO 4/19
* Build a more robust way to grab a subset of a simulation

"""

import matplotlib.pyplot as plt

import skyflux.deprecated.polSims as pol
import skyflux.ant as ant

import numpy as np

def slicer(ant1, ant2, func_show, Qi=None,
    fs=np.arange(125e6, 175e6 + 0.001, 4e6)):
    """
    The second pair of arguments is admittedly
    disappointing, but I could not figure out
    how else to write a script like this.
    """
    # sp doesn't matter since the use of the
    # special_request parameter branches differently
    special = func_show("max_one", Qi=Qi,
        special_request=(ant1, ant2))
    
    b = ant.baselength(ant1, ant2)
    #print(special.shape)
    print("Baseline length [m]:", np.around(b, 4))
    
    center_f = np.average(fs)
    z = pol.fq2z(center_f / 1e9)
    lambda_ = pol.C / center_f
    k_orth = pol.k_perp(z) * b / lambda_
    
    print("Baseline has perpendicular k mode:",
        np.around(k_orth, 4), "[h / Mpc]")
        
    horizonp = pol.horizon_limit(k_orth, z)
    horizonn = -pol.horizon_limit(k_orth, z)
    
    magnitudes = []
    
    # this loop only serves to re-organize the results
    for stokes_param in special:
        magnitudes.append([])
        
        i = len(magnitudes) - 1
        magnitudes[i].append(
            np.fft.fftshift(stokes_param[:, 0]))
        magnitudes[i].append(
            np.fft.fftshift(stokes_param[:, 1]))
        magnitudes[i] = np.array(magnitudes[i])
    
    magnitudes = np.array(magnitudes)
    
    """ OG Bild
    plt.plot(magnitudes[0, 0], magnitudes[0, 1], label="I")
    plt.plot(magnitudes[1, 0], magnitudes[1, 1], label="Q")
    plt.plot(magnitudes[2, 0], magnitudes[2, 1], label="U")
    plt.plot(magnitudes[3, 0], magnitudes[3, 1], label="V")
    """
    
    #plt.plot(magnitudes[0, 0], 10 ** (magnitudes[0, 1] - magnitudes[0, 1]), label="I / I")
    plt.plot(magnitudes[1, 0], 10 ** (magnitudes[1, 1] - magnitudes[0, 1]), label="Q / I")
    plt.plot(magnitudes[2, 0], 10 ** (magnitudes[2, 1] - magnitudes[0, 1]), label="U / I")
    plt.plot(magnitudes[3, 0], 10 ** (magnitudes[3, 1] - magnitudes[0, 1]), label="V / I")
    
    
    """ outdated, but maybe you want to res it to keep debugging
    ### try to find the maximum value associated with
    ### each index
    
    max_counts = [0, 0, 0, 0]
    
    for vi in range(len(special[0, :, yi])):
        vis = [0, 0, 0, 0]
        vis[0] = special[0, vi, yi]
        vis[1] = special[1, vi, yi]
        vis[2] = special[2, vi, yi]
        vis[3] = special[3, vi, yi]
    
        vis_test = np.array([
            np.abs(S) for S in vis
        ])
        
        for si in range(len(vis_test)):
            if vis_test[si] == vis_test.max():
               max_counts[si] += 1 
    
    print(max_counts)
    ###
    """
    
    ### this is pretty bad
    
    ymin = magnitudes[0, 1].min()
    ymin = min(ymin, magnitudes[1, 1].min())
    ymin = min(ymin, magnitudes[2, 1].min())
    ymin = min(ymin, magnitudes[3, 1].min())
    
    ymax = magnitudes[0, 1].max()
    ymax = max(ymin, magnitudes[1, 1].max())
    ymax = max(ymin, magnitudes[2, 1].max())
    ymax = max(ymin, magnitudes[3, 1].max())
    
    ###
    
    """
    plt.vlines(horizonp, ymin, ymax,
        linestyles="dashed", colors='k')
    plt.vlines(horizonn, ymin, ymax,
        linestyles="dashed", colors='k')
    """
    
    plt.xlabel("$k_\parallel$ [$h$ Mpc$^{-1}$]")
    plt.ylabel("log$_{10}$ [K$^2$ ($h^{-1}$ Mpc)^3] ?")
    plt.title(str(ant1) + "-" + str(ant2) + " baseline")
    plt.legend(loc='upper right')
    
    # You should investigate the horizons for
    # any lingering amplitude.
    
    plt.show()
    return special

