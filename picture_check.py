#! /usr/bin/env python

### I should probably kill this script;
### I think it is entirely obsolete in light of pps.py

"power plot sketch"

import matplotlib.pyplot as plt
import numpy as np

import pickle
picture_file = open("67_src_wedge.pickle", "rb")
meta = pickle.load(picture_file)
pic = meta['picture']

for ant1 in pic.keys():
    for ant2 in pic[ant1].keys():
        for nu_idx in range(len(pic[ant1][ant2])):
            next_plot = np.array(pic[ant1][ant2][nu_idx])
            plt.plot(np.abs(next_plot[:, 0]))
            plt.show()

