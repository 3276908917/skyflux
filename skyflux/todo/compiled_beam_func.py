"""
This script's original raison d'etre:

I want to run it once, based on whatever spin 1 harmonics I happen to have
in the repository at the time of running.

Then I want to store this function in a pickle.
When I import skyflux, it should unpickle the function,
    rather than simply recalculating it.
"""

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
beam_origin = os.path.dirname(os.path.abspath(__file__)) + \
              "/HERA_spin1_harmonics.h5"

# My understanding of the documentation is that
# alt, az are both in radians
"""
spline_beam_func = beam_models.model_data_to_spline_beam_func(
    beam_origin,
    # Like in generate_model.py, we have some hard-coded frequencies
    # which we want to re-evaluate in the future.
    beam_frqs
)
"""

"""
nu_axis, tx, ty, kx, ky, E_coeffs, rE_coeffs = beam_models.model_data_to_spline_params(
    beam_origin,
    beam_frqs
)
"""

#! This file is 840 MB. There's no way that's going to fly.
# At the same time, there is no denying that we have shaved a couple of
# seconds off of the import time...
sbfps = utils.load_saves('sbf_params')['arr_0']

spline_beam_func = construct_spline_beam_func(
    sbfps[0], sbfps[1], sbfps[2], sbfps[3],
    sbfps[4], sbfps[5], sbfps[6])

"""

I forgot what I was doing and accidentally wrote the wrong script!!

The whole point of THIS script was to run ONCE. But I have
simply offloaded that part of stokes.py which was supposed to run EVERY
time skyflux gets imported! 

"""
