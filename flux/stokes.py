import os
import numpy as np
from RIMEz import beam_models

from flux import rot

beam_origin = os.path.dirname(os.path.abspath(__file__)) + \
              "/ant.h5"

spline_beam_func = beam_models.model_data_to_spline_beam_func(
    beam_origin,
    # Like in generate_model.py, we have some hard-coded frequencies
    # which we want to re-evaluate in the future.
    np.array([150e6, 151e6, 152e6])
)

S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

def J_matrix(ra, dec, nu=150e6):
    """
    @source : we expect an object of the type
                GLEAM_entry (see parse.py)

    The default argument comes from the beam that I
    had access to when this was written.
    """
    latitude = np.radians(rot.hera_lat)
    # We want to transform the following into a function parameter
    lst = rot.get_lst()
    az, alt = rot.eq_to_topo(ra, dec, latitude, lst, radians=True)

    az = np.array([az])
    alt = np.array([alt])

    #! bad hard-coding
    return spline_beam_func(freq, alt, az)

def A_matrix(ra, dec, nu=150e6):
    J = J_matrix(ra, dec, nu)
    J_outer = np.kron(J, np.conj(J))
    return np.dot(S, np.dot(J_outer, np.linalg.inv(S)))
