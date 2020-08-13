"""
This is a modified form of skyflux's stokes.py and vis.py
It allows for direct input of the J matrix, for debugging purposes
"""

import numpy as np

from skyflux import catalog
from skyflux import rot
from skyflux import ant
from skyflux import stokes

S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

def A_matrix(J):
    J_outer = np.kron(J, np.conj(J))
    return np.dot(S, np.dot(J_outer, np.linalg.inv(S)))

def visibility(J, ant1, ant2, source, ra, dec, nu=151e6, time=None):
    """
    Visibility integrand evaluated for a single source.
    @ant1 and @ant2 are indices of antannae,
        to specify a baseline
    @source is a GLEAM catalog object
        (see catalog.py for specifications)
    @nu : signal frequency [Hz]
    @time : local sidereal time [float, radians]
        default: None corresponds to run-time LST.
    """
    print("time:", time)
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    r = rot.radec2lm(ra, dec, ra0=time)

    phi = ant.phase_factor(ant1, ant2, r, nu)

    malformed_result = np.dot(np.dot(A_matrix(J), s), phi)
    # Hack to get rid of extra array shell surrounding answer
    return malformed_result[:, 0]
