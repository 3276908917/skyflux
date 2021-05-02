"""
TODO 4/19
* Build a more robust way to grab a subset of a simulation

"""


import matplotlib.pyplot as plt

import skyflux.deprecated.polSims as pol
import skyflux.ant as ant

import numpy as np

def slicer(ant1, ant2, func_show, Qi,
    fs=np.arange(50e6, 250e6 + 0.001, 1e6)):
    """
    The latter two arguments are admittedly
    disappointing, but I could not figure out
    how else to write a script like this.
    """
    special = func_show("0E300-387w", 0, Qi=Qi,
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
    
    # you should automatically print out
        # corresponding k_perpendicular
        # baseline length
    
    # Furthermore, you should plot the horizon lines here
    
    # Furthermore, you should investigate the horizons for
    # any lingering amplitude.
    
    plt.show()
    return special

