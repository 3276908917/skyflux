import matplotlib.pyplot as plt

def slicer(ant1, ant2, func_show, Qi):
    """
    The latter two arguments are admittedly
    disappointing, but I could not figure out
    how else to write a script like this.
    """
    special = func_show("E0-387w", 0, Qi=Qi,
        special_request=(ant1, ant2))
    
    print(special.shape)    
        
    plt.plot(special[0, :, 0], special[0, :, 1], label="I")
    plt.plot(special[1, :, 0], special[1, :, 1], label="Q")
    plt.plot(special[2, :, 0], special[2, :, 1], label="U")
    plt.plot(special[3, :, 0], special[3, :, 1], label="V")
    
    plt.xlabel("$k_\parallel$ [$h$ Mpc$^{-1}$]")
    plt.ylabel("log$_{10}$ [K$^2$ ($h^{-1}$ Mpc)^3] ?")
    plt.title(str(ant1) + "-" + str(ant2) + " baseline")
    
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

