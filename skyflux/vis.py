"""
Utilities for evaluating the visibility sky integral.
"""

import numpy as np

from skyflux import rot
from skyflux import ant
from skyflux import stokes
from skyflux import catalog

#! the default arguments for 'nu' are inconsistent across functions...
def visibility(ant1, ant2, source, nu=151e6, time=None):
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
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    A = stokes.create_A(ra=ra, dec=dec, lst=time, radians=True)
    
    r = rot.radec2lm(ra, dec, ra0=time)
    phi = ant.phase_factor(ant1, ant2, r, nu)

    malformed_result = np.dot(np.dot(A, s), phi)
    # Hack to get rid of extra array shell surrounding answer
    return malformed_result[:, 0]


def sources_over_time(ant1, ant2, list_sources=None,
                        start=0, end=2/3*np.pi, interval=np.pi/72, nu=151e6):
    """
    Return an array containing the visibilities at different points of time.
    @ant1 and @ant2 are indices of antennae, to specify a baseline
    @list_sources is an array of GLEAM catalog objects (see catalog.py for specifications)
        default value translates to the entire downloaded segment of the catalog
    @start: starting LST of integration [float, radians]
        default: 0 hours (cold patch)
    @end: terminal LST of integration [float, radians]
        default: 8 hours = 2/3 pi (cold patch)
    @interval: integration window width [float, radians]
        default: 10 minutes = np.pi / 72
    @nu frequency in Hertz
    """
    if list_sources is None
        # make a copy to ensure write safety
        list_sources = catalog.obj_catalog.copy()

    # If the user has entered a single source directly,
    # we can automatically standardize the formatting
    # by placing it by itself in a list
    if type(list_sources) != list:
        list_sources = [list_sources]
    
    list_visibilities = []
    lst = start
    while lst <= end:
        next_vista = np.array([0j, 0j, 0j, 0j])
        for source in list_sources:
            next_vista += visibility(ant1, ant2, source, nu, time)

        list_visibilities.append(next_vista)
        lst += interval
    # perhaps not necessary. Better safe than sorry:
    return np.array(list_visibilities)
