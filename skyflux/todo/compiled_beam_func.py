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

nu_axis, tx, ty, kx, ky, E_coeffs, rE_coeffs = beam_models.model_data_to_spline_params(
    beam_origin,
    beam_frqs
)

sbf_storage = np.array([nu_axis, tx, ty, kx, ky, E_coeffs, rE_coeffs])

import numpy as np
import os
from RIMEz import beam_models
from skyflux import utils

# disgusting hack
MACRO_EPSILON = 0.001

#! I guess, to stay on the safe side, we will want to pickle
    # the beam_frequencies separately from the spline_beam_func. 
beam_frqs = np.arange(50e6, 250e6 + MACRO_EPSILON, 1e6)

# It is imperative that the name here line up with that used in
# generate_model.py (a file NOT included in the installation, but
# included in the Git repository!)
sbfps_origin = os.path.dirname(os.path.abspath(__file__)) + \
              "/sbf_params.npz"



"""

I forgot what I was doing and accidentally wrote the wrong script!!

The whole point of THIS script was to run ONCE. But I have
simply offloaded that part of stokes.py which was supposed to run EVERY
time skyflux gets imported! 

"""
