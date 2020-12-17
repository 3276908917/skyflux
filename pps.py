"power plot sketch"

import skyflux as sf

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

import pickle

picture_file = open("picture_dict.pickle", "rb")

meta = pickle.load(picture_file)

frq = meta['frequencies']
ap = meta['ant_pos']
pic = meta['picture']

import copy
avg_pic = pic.copy() # terrible on RAM, but automatically provides a skeleton

for ant1 in pic.keys():
    for ant2 in pic[ant1].keys():
        """
        Since we only used one LST,
        we do not need to do any averaging
            (averaging is supposed to happen over LSTs, not frequency)

        If you wanted to stick with your current data,
            you can take the norm squared of those single LSTs,
            but in the real world that would bring on heaps of noise.
        """

        # this is a proportionality.
        # The real deal uses the power equation 6 from Nunhokee et al.
        power = np.vdot(pic[ant1][ant2][nu], pic[ant1][ant2][nu + 1])
        avg_pic[ant1][ant2] = power

"""
Effective baseline is 200 MHz? Since I am running from 50-250 MHz.
"""


"""
These functions are obsolete:
    polarizedSims has a COSMO_constants script which can handle these
"""

# It's not lmn space after all... it is lm eta space

def orthk(D, lambda_, ant1, ant2):
    """
    ?? D: transverse comoving distance
    lambda_: observing wavelength
    """
    # I am guessing I get lambda_ by using c = lambda * f where
    # f is the beam frequency
    b_norm = np.linalg.norm(ant.get_baseline(ant1, ant2))
    return 2 * np.pi * b_norm / lambda_ / D

# magic numbers to directly transcribe paper
H_0 = 100 # h km / s what is h here?
Omega_M = 0.27
Omega_k = 0
Omega_Lambda = 0.73

#? 21-cm line rest frequency
f21 = 1420405721.7667
def parak(z, ant1, ant2):
    """
    z: redshift
    """
    # unfortunately, this gives us baseline in UVW, not lm\eta coordinates
    b = ant.get_baseline(ant1, ant2)
    eta = something(b)
    
    Z = z + 1
    E = (Omega_M * Z ** 3 + Omega_k * Z ** 2 + Omega_Lambda) ** 0.5
    return eta * np.pi * 2 * f21 * H_0 * E / c / Z ** 2
