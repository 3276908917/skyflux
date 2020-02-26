# All credit to https://pyuvdata.readthedocs.io/en/latest/tutorial.html

import os
from pyuvdata import UVBeam
from pyuvdata.data import DATA_PATH
import numpy as np
import matplotlib.pyplot as plt
beam = UVBeam()

beam_filenames = [os.path.join(DATA_PATH, f) for f in ['/home/lfinkbeiner/Documents/HERA/HERA_4.9m_E-pattern_151MHz.txt']]

# Most importantly, the hard-code values below are entirely taken from the demo,
    # and do not necessarily reflect your case at all
# It seems a little ridiculous to have to specify frequencies both here and in the file names...
# I do not understand feed polarization. And what is the point of feed version?
beam.read_cst_beam(beam_filenames, beam_type='efield', frequency=[151e6],
    feed_pol='x', rotate_pol=True, telescope_name='HERA',
    feed_name='PAPER_dipole', feed_version='0.1',
    model_name='E-field pattern - Rigging height 4.9m',
    model_version='1.0')
print("Beam type", beam.beam_type)
print("Pixel coordinate system", beam.pixel_coordinate_system)
# I do not know what this means
print("Data normalization", beam.data_normalization)

print("Number of beam polarizations", beam.Npols)
print("Polarization type", beam.polarization_array)
# this should be one, in my case
print("Number of input frequencies", beam.Nfreqs)
# I am not gathering anything helpful from this either
print("Array data shape", beam.data_array.shape)

# "plot zenith angle cut through beam" I am not sure how relevant this is to my current effort
plt.plot(beam.axis2_array, np.real(beam.data_array[0, 0, 0, 0, :, 0]), label='real') # fancy syntax, what is going on here?
plt.plot(beam.axis2_array, np.imag(beam.data_array[0, 0, 0, 0, :, 0]), label='imaginary')
plt.xscale('log')
plt.xlabel('Zenith Angle (radians)')
plt.ylabel('Power')
plt.legend(bbox_to_anchor=(1,1))

# I do not understand how power can have an imaginary component.
plt.show()

print(beam.data_array[0, 0, 0, 0, :, 0])
# hey, I didn't mean anything
