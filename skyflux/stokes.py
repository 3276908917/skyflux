import os
import numpy as np
from RIMEz import beam_models

from skyflux import rot

beam_origin = os.path.dirname(os.path.abspath(__file__)) + "/ant.h5"

# My understanding of the documentation is that
# alt, az are both in radians
spline_beam_func = beam_models.model_data_to_spline_beam_func(
    beam_origin,
    # Like in generate_model.py, we have some hard-coded frequencies
    # which we want to re-evaluate in the future.
    np.array([150e6, 151e6, 152e6])
)

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

def create_J(ra=None, dec=None, az=None, alt=None,
             lat=None, lst=None, nu=150e6):
    """
    Return the Jones matrix J.
    @ra: right ascension of the source, in radians.
    @dec: right ascension of the source, in radians.
    @lat: latitude of point of observation, in radians
        default: HERA array
    @lst: local sidereal time, in radians
        default: time of execution
    
    @nu: frequency of interest, in Hz.
        default: 150 MHz
    The default argument comes from the beam that I
    had access to when this was written.
    """
    # This section handles the many possible bad inputs
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
        az, alt = rot.eq_to_topo(ra, dec, lat=lat, lst=lst, radians=True)

    # This section handles different input objects
    if type(az) == list:
        az = np.array(az)
    if type(alt) == list:
        alt = np.array(alt)

    if type(az) != np.ndarray and type(alt) != np.ndarray
        az = np.array([az])
        alt = np.array([alt])
    elif type(az) == np.ndarray and type(az) == np.ndarray:
        return format_J(spline_beam_func(nu, alt, az))
    
    raise TypeError('shapes of inputs must be the same (got one list and one scalar).')

def create_A(ra=None, dec=None, az=None, alt=None, J=None,
             lat=None, lst=None, nu=150e6):
    """
    Return the Mueller matrix A.
    @ra: right ascension of the source, in radians.
    @dec: declination of the source, in radians.
    @lat: latitude of point of observation, in radians
        default: HERA array
    @lst: local sidereal time, in radians
        default: time of execution
        
    @nu: frequency of interest, in Hz.
        default: 150 MHz
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
        J = J_matrix(ra, dec, az, alt, lat, lst, nu)

    # It would be cool to switch functionality if we were given an array of J's
    # in which case, we should individually calculate the A matrix for each one, and return that.
    
    J_outer = np.kron(J, np.conj(J))
    return np.dot(S, np.dot(J_outer, np.linalg.inv(S)))
