import os
import numpy as np
from RIMEz import beam_models

from flux import rot

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
    for J in new_J:
        J[0, 0], J[1, 0], J[0, 1], J[1, 1] = \
        J[1, 0], J[0, 0], J[1, 1], J[0, 1]
    return new_J

# This is a constant change of basis matrix
# for getting stokes parameters with a Jones matrix.
S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

def J_matrix(ra, dec, lat=None, lst=None, nu=150e6):
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
    az, alt = rot.eq_to_topo(ra, dec, lat=lat, lst=lst, radians=True)
    az = np.array([az])
    alt = np.array([alt])

    return format_J(spline_beam_func(nu, alt, az))

def A_matrix(ra, dec, lat=None, lst=None, nu=150e6):
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
    J = J_matrix(ra, dec, lat, lst, nu)
    J_outer = np.kron(J, np.conj(J))
    return np.dot(S, np.dot(J_outer, np.linalg.inv(S)))
