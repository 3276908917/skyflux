"""
This is currently just a placeholder script.

What do I want it to do?

I want to run it once, based on whatever GLEAM catalog files I happen
to have in the repository at the time of running.

It will generate a numpy array of GLEAM objects,
than I will use numpy savez to get the array into a disk file.
"""

import numpy as np

from skyflux import catalog

import os
data_prefix = os.path.dirname(os.path.abspath(__file__)) + "/"

try:
    f = open(data_prefix + "gleam_with_alpha.txt", "r")
    obj_catalog = []
    # For each line in f, the delimiter is |
    for line in f:
        obj_catalog.append(catalog.GLEAM_entry(line[1:]))
    f.close()
    np.save(data_prefix + '../catalog.npy', obj_catalog, allow_pickle=True)
except FileNotFoundError:
    print("Failure to load gleam catalog.")
