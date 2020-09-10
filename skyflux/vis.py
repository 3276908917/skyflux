"""
Utilities for evaluating the visibility sky integral.
"""

import numpy as np

from skyflux import rot
from skyflux import ant
from skyflux import stokes
from skyflux import catalog

#! the default arguments for 'nu' are inconsistent across functions...
def visibility(ant1, ant2, source, nu=151, time=None):
    """
    Visibility integrand evaluated for a single source.
    @ant1 and @ant2 are indices of antannae,
        to specify a baseline
    @source is a GLEAM catalog object
        (see catalog.py for specifications)
    @nu : signal frequency [MHz]
    @time : local sidereal time [float, radians]
        default: None corresponds to run-time LST.
    """
    # make a copy to ensure write safety
    ef = catalog.expected_frequencies.copy()
    
    # keep in mind that ef is in ascending order
    if nu in ef: # we have an observed value for this frequency
        I = source.flux_by_frq[nu]
    elif nu < ef[0] or nu > ef[len(ef)-1]: # "danger of extrapolation" --Aaron Simon
        raise NotImplementedError("That frequency would have to be extrapolated.")
    else: # we linearly interpolate between the two neighboring frequencies
        # I imagine there are several ways to refine this approach,
            # (for example, include the frequency before last and the frequency after next, too)
        # and linear interpolation may not even be a reasonable model for this problem,
        # so take this suite as a placeholder
        i = 0

        # Use power law formula here; use the spectral index
        # s1/s2 = nu1/nu2 * s2
        while expected_frequencies[i] > nu:
            i += 1
        nu_a = expected_frequencies[i - 1]
        nu_b = expected_frequencies[i]
        span = nu_b - nu_a
        interp_a = (nu_b - nu) * source.flux_by_frq[nu_a]
        interp_b = (nu - nu_a) * source.flux_by_frq[nu_b]
        I = (interp_a + interp_b) / span
        
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
    @ant1 and @ant2 are indices of antennae, to specify a baseline.
    @list_sources is an array of GLEAM catalog objects (see catalog.py for specifications)
        default value translates to the entire downloaded segment of the catalog.
    @start: starting LST of integration [float, radians]
        default: 0 hours (cold patch)
    @end: terminal LST of integration [float, radians]
        default: 8 hours = 2/3 pi (cold patch)
    @interval: integration window width [float, radians]
        default: 10 minutes = np.pi / 72
    @nu frequency in Hertz
    """
    if list_sources is None:
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
            next_vista += visibility(ant1, ant2, source, nu=nu, time=lst)

        list_visibilities.append(np.array([lst, next_vista]))
        lst += interval
    # perhaps not necessary. Better safe than sorry:
    return np.array(list_visibilities)

def sources_over_frequency(ant1, ant2, list_sources=None,
                        start=76e6, end=227e6, interval=1e6, time=None):
    """
    Return an array containing the visibilities at different values for the frequency
        of incident light.
    @ant1 and @ant2 are indices of antennae, to specify a baseline.
    @list_sources is an array of GLEAM catalog objects (see catalog.py for specifications)
        default value translates to the entire downloaded segment of the catalog.
    @start: starting frequency of integration 
        default: 76 MHz (minimum frequency described by the GLEAM catalog)
    @end: terminal frequency of integration
        default: 227 MHz
    @interval: integration window width
        default: 1 MHz (arbitrary)
    @time fixed LST at which to evaluate all visibilities
        default: None translates to the LST as of the initial call of this function
    """
    if list_sources is None:
        # make a copy to ensure write safety
        list_sources = catalog.obj_catalog.copy()
    if time is None:
        time = rot.get_lst(radians=True)

    # If the user has entered a single source directly,
    # we can automatically standardize the formatting
    # by placing it by itself in a list
    if type(list_sources) != list:
        list_sources = [list_sources]
    
    list_visibilities = []
    lfreq = start
    while lfreq <= end:
        next_vista = np.array([0j, 0j, 0j, 0j])
        for source in list_sources:
            next_vista += visibility(ant1, ant2, source, nu=lfreq, time=time)

        list_visibilities.append(np.array([lfreq, next_vista]))
        lfreq += interval
    # perhaps not necessary. Better safe than sorry:
    return np.array(list_visibilities)
