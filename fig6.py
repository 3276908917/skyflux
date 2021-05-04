"""
TODO 4/19
* Build a more robust way to grab a subset of a simulation

"""

import matplotlib.pyplot as plt

import skyflux.deprecated.polSims as pol
import skyflux.ant as ant

import numpy as np

def slicer(ant1, ant2, func_show, Qi=None,
    fs=np.arange(50e6, 250e6 + 0.001, 4e6)):
    """
    The second pair of arguments is admittedly
    disappointing, but I could not figure out
    how else to write a script like this.
    """
    #!!! Isn't it kind of weird that we're using
    # sp=0 here?
    special = func_show("0E300-387w", sp=0, Qi=Qi,
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
    
    special = np.fft.fftshift(special)
    
    yi = 0
    xi = 1
    
    plt.plot(special[0, :, xi], special[0, :, yi], label="I")
    plt.plot(special[1, :, xi], special[1, :, yi], label="Q")
    plt.plot(special[2, :, xi], special[2, :, yi], label="U")
    plt.plot(special[3, :, xi], special[3, :, yi], label="V")
    
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
    
    ### this is pretty bad
    
    ymin = special[0, :, yi].min()
    ymin = min(ymin, special[1, :, yi].min())
    ymin = min(ymin, special[2, :, yi].min())
    ymin = min(ymin, special[3, :, yi].min())
    
    ymax = special[0, :, yi].max()
    ymax = max(ymin, special[1, :, yi].max())
    ymax = max(ymin, special[2, :, yi].max())
    ymax = max(ymin, special[3, :, yi].max())
    
    ###
    
    plt.vlines(horizonp, ymin, ymax,
        linestyles="dashed", colors='k')
    plt.vlines(horizonn, ymin, ymax,
        linestyles="dashed", colors='k')
    
    plt.xlabel("$k_\parallel$ [$h$ Mpc$^{-1}$]")
    plt.ylabel("log$_{10}$ [K$^2$ ($h^{-1}$ Mpc)^3] ?")
    plt.title(str(ant1) + "-" + str(ant2) + " baseline")
    plt.legend(loc='upper right')
    
    # You should investigate the horizons for
    # any lingering amplitude.
    
    plt.show()
    return special

