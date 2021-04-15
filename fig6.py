import matplotlib.pyplot as plt

import skyflux.deprecated.polSims as pol
import skyflux.ant as ant

import numpy as np

def slicer(ant1, ant2, func_show, Qi,
    fs=np.arange(125e6, 175e6 + 0.001, 1e6)):
    """
    The latter two arguments are admittedly
    disappointing, but I could not figure out
    how else to write a script like this.
    """
    special = func_show("E0-387w", 0, Qi=Qi,
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
        
    plt.plot(special[0, :, 0], special[0, :, 1], label="I")
    plt.plot(special[1, :, 0], special[1, :, 1], label="Q")
    plt.plot(special[2, :, 0], special[2, :, 1], label="U")
    plt.plot(special[3, :, 0], special[3, :, 1], label="V")
    
    ### this is pretty bad
    
    ymin = special[0, :, 1].min()
    ymin = min(ymin, special[1, :, 1].min())
    ymin = min(ymin, special[2, :, 1].min())
    ymin = min(ymin, special[3, :, 1].min())
    
    ymax = special[0, :, 1].max()
    ymax = max(ymin, special[1, :, 1].max())
    ymax = max(ymin, special[2, :, 1].max())
    ymax = max(ymin, special[3, :, 1].max())
    
    ###
    
    plt.vlines(horizonp, ymin, ymax,
        linestyles="dashed", colors='k')
    plt.vlines(horizonn, ymin, ymax,
        linestyles="dashed", colors='k')
    
    plt.xlabel("$k_\parallel$ [$h$ Mpc$^{-1}$]")
    plt.ylabel("log$_{10}$ [K$^2$ ($h^{-1}$ Mpc)^3] ?")
    plt.title(str(ant1) + "-" + str(ant2) + " baseline")
    plt.legend(loc='upper right')
    
    # you should automatically print out
        # corresponding k_perpendicular
        # baseline length
    
    # Furthermore, you should plot the horizon lines here
    
    # Furthermore, you should investigate the horizons for
    # any lingering amplitude.
    
    plt.show()
    return special
    
"""
    for i in range(len(horizonx)):
        k_orthogonal = horizonx[i]
        zi = int(i * len(fs) / len(horizonx)) # terrible
        zloc = pol.fq2z(fs[zi] / 1e9)
        # after fixing this: pick a bin of k_perpendicular and do a cut
        # it's a 1D plot, include all four Stokes parameters
        # as different lines
        baselength = k_orthogonal / k_starter
        print(baselength, "baselength [m]")
        tau = baselength / 2.9979e8 # geometric delay, Parsons eq. 3
        print(tau * 1e9, "delay [ns]")
        horizonyp.append(pol.k_parallel(tau, zloc))
        horizonym.append(pol.k_parallel(-tau, zloc))
    """

