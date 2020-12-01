import numpy as np
import healpy as hp

from skyflux import rot

from skyflux.todo.compiled_beam_func import spline_beam_func
from skyflux.todo.compiled_beam_func import beam_frqs

# disgusting hack
MACRO_EPSILON = 0.001

def format_J(J_RIMEz):
    """
    The default shape of a single RIMEz J matrix is:
        [[xy, yy], [xx, yx]]
    The true convention is:
        [[xx, xy], [yx, yy]]
    As with np.radians, one must verify that the input is of
    the incorrect format. For example, if one called this on
    a hard-coded J matrix, one would merely confuse the order.
    """
    new_J = np.copy(J_RIMEz)
    for i in range(len(new_J)):
        new_J[i, 0, 0], new_J[i, 1, 0], new_J[i, 0, 1], new_J[i, 1, 1] = \
        new_J[i, 0, 1], new_J[i, 0, 0], new_J[i, 1, 1], new_J[i, 1, 0]
    return new_J

# This is a constant change of basis matrix
# for getting stokes parameters with a Jones matrix.
S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])
Si = np.linalg.inv(S)

def create_J(ra=None, dec=None, az=None, alt=None,
             lat=None, lst=None, nu=151e6, radians=False):
    """
    Return the Jones matrix J.
    @ra: right ascension of the source, in radians.
    @dec: right ascension of the source, in radians.
    @lat: latitude of point of observation, in radians
        default: HERA array
    @lst: local sidereal time, in radians
        default: time of execution
    
    @nu: frequency of interest, in Hz.
        default: 151 MHz
    The default argument comes from the beam that I
    had access to when this was written.
    """
    # This section handles the many possible bad combinations of inputs
    if type(nu) == list:
        nu = np.array(list)
    
    if type(nu) == np.ndarray:
        for frequency in nu:
            if frequency not in beam_frqs:
                raise NotImplementedError("No routine for interpolating between beam frequencies.")
    elif nu not in beam_frqs:
        raise NotImplementedError("No routine for interpolating between beam frequencies.")
    
    if ra is not None:
        if dec is None:
            raise TypeError('ra was provided without accompanying dec.')
        if az is not None:
            raise TypeError('create_J accepts only one coordinate system; both ra and az were given.')
        if alt is not None:
            raise TypeError('create_J accepts only one coordinate system; both ra and alt were given.')
    if dec is not None and ra is None:
        raise TypeError('dec was provided without accompanying ra.')   
    if az is not None and alt is None:
        raise TypeError('az was provided without accompanying alt.')
    if alt is not None and az is None:
        raise TypeError('alt was provided without accompanying az.')

    # Now we ensure we are in topocentric coordinates
    if az is None and alt is None:
        az, alt = rot.eq_to_topo(ra, dec, lat=lat, lst=lst, radians=radians)

    if not radians:
        az = np.radians(az)
        alt = np.radians(alt)

    # This section handles different input objects
    if type(az) == list:
        az = np.array(az)
    if type(alt) == list:
        alt = np.array(alt)

    if type(az) != np.ndarray and type(alt) != np.ndarray:
        az = np.array([az])
        alt = np.array([alt])

    # one of the inputs did not get converted
    if type(az) != np.ndarray or type(alt) != np.ndarray:
        if type(az) != np.ndarray:
            az = az * np.ones(len(alt))
        elif type(alt) != np.ndarray:
            alt = alt * np.ones(len(alt))
        #raise TypeError('shapes of inputs must be the same (got one list and one scalar).')    
        
    return format_J(spline_beam_func(nu, alt, az))

def create_J_sky(nside, nu=151e6):
    """
    Create a deck of J matrices to describe the whole sky
    @nside : integer, power of 2
        greater nside yields a more finely-grained deck
    @nu : Hz
        frequency of beam for which we are generating a whole-sky deck
    Use create_A_space if a Mueller sky is desired as the output.
    """
    theta, phi = hp.pix2ang(nside, np.arange(12 * nside ** 2))
    az = phi
    alt = np.pi / 2 - theta
    J_raw = spline_beam_func(nu, alt, az)
    return format_J(J_raw)

create_A_sky = lambda nside, nu=151e6: \
    np.array([create_A(J=Ji) for Ji in create_J_space(nside, nu)])

def create_A(ra=None, dec=None, az=None, alt=None, J=None,
             lat=None, lst=None, nu=151e6, radians=False):
    """
    Return the Mueller matrix A.
    @ra: right ascension of the source, in radians.
    @dec: declination of the source, in radians.
    @lat: latitude of point of observation, in radians
        default: HERA array
    @lst: local sidereal time, in radians
        default: time of execution
        
    @nu: frequency of interest, in Hz.
        default: 151 MHz
    The default argument comes from the beam that I
    had access to when this was written.
    """
    if J is not None:
        if ra is not None:
            raise TypeError('create_A accepts only one input representation; both J and ra were given.')
        if dec is not None:
            raise TypeError('create_A accepts only one input representation; both J and dec were given.')
        if az is not None:
            raise TypeError('create_A accepts only one input representation; both J and az were given.')
        if alt is not None:
            raise TypeError('create_A accepts only one input representation; both J and alt were given.')
    else:
        J = create_J(ra, dec, az, alt, lat, lst, nu, radians)

    # It would be cool to switch functionality if we were given an array of J's
    # in which case, we should individually calculate the A matrix for each one, and return that.
    
    J_outer = np.kron(J, np.conj(J))
    return np.dot(Si, np.dot(J_outer, S))
