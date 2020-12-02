"""
We want to run this script once, based on whatever spin 1 harmonics
    I happen to have in the repository at the time of running.

Then I want to store the spline beam func parameters in a numpy
    array to file.
When I import skyflux, it should load the array and reconstruct
    the spline beam function from these parameters.

Originally, I wanted to store the function itself somewhere, but this
    is the closest thing that I know how to do
(dill was certainly unable to pickle the spline beam function).
"""

print("\nStarting script to generate spline beam function " + \
      "parameters according to the spin-1 harmonics in this directory.")

import numpy as np
import os
from RIMEz import beam_models

here = os.path.dirname(os.path.abspath(__file__)) + "/"

# It is imperative that the names here line up with those used in
# generate_model.py (a file NOT included in the installation, but
# included in the Git repository!)
beam_origin = here + "HERA_spin1_harmonics.h5"
frqs_origin = "HERA_beam_frqs.npy"

beam_frqs = np.load(here + frqs_origin)

print("\nAll header variables processed. Computing SBF params...")

nu_axis, tx, ty, kx, ky, E_coeffs, rE_coeffs = \
    beam_models.model_data_to_spline_params(
        beam_origin,
        beam_frqs
    )

print("\nParams computed. Storing full results...")

sbf_storage = np.array([nu_axis, tx, ty, kx, ky, E_coeffs, rE_coeffs])

np.save(here + "sbf_params.npz", sbf_storage, allow_pickle=True)
