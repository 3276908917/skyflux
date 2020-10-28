"""
Utilities for evaluating the visibility sky integral.
"""

import numpy as np

from skyflux import rot
from skyflux import ant
from skyflux import stokes
from skyflux import catalog

# Spectral index, coefficient of proportionality, empirically tuned
cop = lambda source, nu : source.flux_by_frq[nu / 1e6] / nu ** source.alpha

def get_I(source, nu=151e6):
    # make a copy to ensure write safety
    ef = catalog.expected_frequencies.copy()

    index_nu = nu / 1e6
    
    # keep in mind that ef is in ascending order
    #if index_nu in ef: # we have an observed value for this frequency
    #    return source.flux_by_frq[index_nu]

    #elif index_nu < ef[0] or index_nu > ef[len(ef)-1]:
        # "beware the danger of extrapolation" --Aaron Simon
        #raise NotImplementedError("That frequency would have to be extrapolated.")

    if index_nu < ef[0] or index_nu > ef[len(ef)-1]:
        raise UserWarning("Using an extrapolated frequency: " + str(nu))

    # we use the power law spectral index formulation here
        # conceptually, I may have made an error.
        # I think I want flux, but I am here calculating flux per unit frequency
        # See "single_source_over_nu.ipynb" for an example plot of how alpha varies...
    i = 0

    # Use power law formula here; use the spectral index
    # s1/s2 = nu1/nu2 * s2
    while ef[i] > index_nu:
        i += 1

    # You should double-down on interpolation;
    # give up on trying to use any real values for flux
    # we want a smooth spectrum, and these real references
    # create discontinuities.

    nu_a = ef[i - 1] * 1e6
    nu_b = ef[i] * 1e6
    span = nu_b - nu_a

    cop_a = cop(source, nu_a)
    cop_b = cop(source, nu_b)

    interp_a = (nu_b - nu) * cop_a
    interp_b = (nu - nu_a) * cop_b
    this_cop = (interp_a + interp_b) / span

    return this_cop * nu ** source.alpha

def visibility(ant1, ant2, source, nu=151e6, time=None):
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
    I = get_I(source, nu)
        
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    
    A = stokes.create_A(ra=ra, dec=dec, lst=time, nu=nu, radians=True)
    r = rot.radec2lm(ra, dec, ra0=time)
    
    phi = ant.phase_factor(ant1, ant2, r, nu)

    malformed_result = np.dot(np.dot(A, s), phi)
    # Hack to get rid of extra array shell surrounding answer
    return malformed_result[:, 0]

# Incoming function, intended to replace sources_over_time
def new_sources_over_time(ant1, ant2, list_sources=None,
                        start=0, end=2/3*np.pi, interval=np.pi/72, nu=151e6, interpolator=None):
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
    if interpolator is None:
        raise NotImplementedError("We are currently hard-coding interpolators.")
    
    if list_sources is None:
        # make a copy to ensure write safety
        list_sources = catalog.obj_catalog.copy()

    # If the user has entered a single source directly,
    # we can automatically standardize the formatting
    # by placing it by itself in a list
    if type(list_sources) != list:
        list_sources = [list_sources]

    list_lst = []
    lst = start
    while lst <= end:
        list_lst.append(lst)
        lst += interval
    list_lst = np.array(list_lst)

    list_visibilities = np.zeros((len(list_lst), 4), dtype=np.complex128)

    for source in list_sources:
        # establish values common to all visibility calculations
        I = get_I(source, nu)
        s = np.array([complex(I), 0, 0, 0])
        ra = np.radians(source.ra_angle)
        dec = np.radians(source.dec_angle)

        for k in range(len(list_lst)):
            curr_lst = list_lst[k]

            # no input for latitude. Code is no longer general :(
            # no input for radians either. We have to rely on the user
                # to provide a radian-based interpolator.
            az, alt = rot.eq_to_topo(ra, dec, lst=curr_lst, radians=True)
            A = interpolator(az, alt)

            r = rot.radec2lm(ra, dec, ra0=curr_lst)
            phi = ant.phase_factor(ant1, ant2, r, nu)

            malformed_result = np.dot(np.dot(A, s), phi)

            list_visibilities[k] += malformed_result
            
    # perhaps not necessary. Better safe than sorry:
    return list_lst, np.array(list_visibilities)

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
    @start: starting frequency of integration, inclusive
        default: 76 MHz (minimum frequency described by the GLEAM catalog)
    @end: terminal frequency of integration, inclusive
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
