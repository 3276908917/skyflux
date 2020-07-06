"""
Utilities for evaluating the visibility sky integral.
"""

import numpy as np

from flux import rot
from flux import ant
from flux import stokes
from flux import catalog

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
    r = rot.raddec2lm(ra, dec, ra0=time)

    phi = ant.phase_factor(ant1, ant2, r, nu)

    malformed_result = np.dot(np.dot(stokes.A_matrix(ra, dec, nu), s), phi)
    # There is something fishy about this array hack.
    # Its necessity bespeaks ill omen
    return malformed_result[:, 0]

def visibility_integrand(ant1, ant2, nu=151e6, time=None):
    """
    Return visibility integrand evaluated for all sources
    in catalog's parsed GLEAM array, obj_catalog.
    @ant1 and @ant2 are indices of antannae,
        to specify a baseline
    @nu frequency in Hertz
    @time : local sidereal time [float, radians]
        default: None corresponds to run-time LST.
    """
    total = np.array([0j, 0j, 0j, 0j]) # 4 x 1. Visibility has a phase,
    for source in catalog.obj_catalog:
        total += visibility(ant1, ant2, source, nu, time)
    return total

def sky_over_time(ant1, ant2, start, end, interval, nu=151e6):
    """
    Return an array containing multiple visibility_integrand results.
    @ant1 and @ant2 are indices of antannae,
        to specify a baseline
    @start: starting LST of integration [float, radians]
    @end: terminal LST of integration [float, radians]
    @interval: integration window width [float, radians]
    @nu frequency in Hertz
    
    The cold patch has start = 0 hours, end = 8 hours,
        and we set out to measure at ten minute intervals.
    (ant1, ant2, 0, 2 / 3 * np.pi, np.pi / 72, nu)
    to-do: error-checking on inputs
    """
    list_visibilities = []
    lst = start
    while lst <= end:
        list_visibilities.append(visibility_integrand(ant1, ant2, nu, lst))
        lst += interval
    # perhaps not necessary. Better safe than sorry:
    return np.array(list_visibilities)

#! Maybe we can merge this function and the previous one somehow,
    # we are re-using virtually all of the code
# chiefly for debugging and discrete analysis purposes
def source_over_time(ant1, ant2, source, start, end, interval, nu=151e6):
    """
    Return an array containing multiple visibility results.
    @ant1 and @ant2 are indices of antannae,
        to specify a baseline
    @source is a GLEAM catalog object
        (see catalog.py for specifications)
    @start: starting LST of integration [float, radians]
    @end: terminal LST of integration [float, radians]
    @interval: integration window width [float, radians]
    @nu frequency in Hertz
    
    The cold patch has start = 0 hours, end = 8 hours,
        and we set out to measure at ten minute intervals.
    (ant1, ant2, 0, 2 / 3 * np.pi, np.pi / 72, nu)
    to-do: error-checking on inputs
    """
    list_visibilities = []
    lst = start
    while lst <= end:
        list_visibilities.append(visibility(ant1, ant2, source, nu, lst))
        lst += interval
    # perhaps not necessary. Better safe than sorry:
    return np.array(list_visibilities)
