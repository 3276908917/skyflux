import pickle
import os

import numpy as np
from RIMEz import beam_models
from spin1_beam_model import cst_processing, jones_matrix_field

data_prefix = os.path.dirname(os.path.abspath(__file__)) + "/"
beam_origin = data_prefix + "HERA_4.9m_E-pattern_151MHz.txt"
beam_destination = data_prefix + "ant1_s2.h5"

processor = cst_processing.CSTDataProcessor(
    [beam_origin, beam_origin, beam_origin],
    # hard coding to demonstrate functionality elsewhere.
    # We will fix it eventually.
    np.array([150e6, 151e6, 152e6]),
    1, 1e-4
)

# If any of the three data files (J matrix, source catalog, antannae positions)
# fails to load, we give a more specific error message for last two functions.
full_load = True

S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

c = 299792458 # m / s

try:
    # Python 2 deserves to die
    J = np.load(data_prefix + "J.npy", fix_imports=False)
    print("Do not forget that our current approach to J" + \
          " pre-generation is imprecise and mostly wrong.")

    def A(source_idx):
        """ Return the Müller matrix for an associated Jones matrix"""
        this_J = J[source_idx]
        J_outer = np.kron(J, np.conj(J))
        return np.dot(np.dot(np.linalg.inv(S), J_outer), S)
except FileNotFoundError:
    print("Failure to load pre-generated J matrix.")
    full_load = False
