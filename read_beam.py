# All credit to https://pyuvdata.readthedocs.io/en/latest/tutorial.html

import os
from pyuvdata import UVBeam
from pyuvdata.data import DATA_PATH
import numpy as np
import matplotlib.pyplot as plt
beam = UVBeam()

beam_filenames = [os.path.join(DATA_PATH, f) for f in ['HERA_4.9m_E-pattern_151MHz.txt']]

beam.read_cst_beam(beam_filenames, beam_type='efield')
print("Beam type", beam.beam_type)
print("Pixel coordinate system", beam.pixel_coordinate_system)
# I do not know what this means
print("Data normalization", data_normalization)

print("Number of beam polarizations", beam.Npols)
print("Polarization type", beam.polarization_array)
# this should be one, in my case
print("Number of input frequencies", beam.Nfreqs)
# I am not gathering anything helpful from this either
print("Array data shape", beam.data_array.shape)

# "plot zenith angle cut through beam" I am not sure how relevant this is to my current effort
plt.plot(beam.axis2_array, beam.data_array[0, 0, 0, 0, :, 0]) # fancy syntax, what is going on here?
plt.xscale('log')
plt.xlabel('Zenith Angle (radians)')
plt.ylabel('Power')
plt.show()

