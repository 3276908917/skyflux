#! /usr/bin/env python

"power plot sketch"

import matplotlib.pyplot as plt
import numpy as np

import pickle
picture_file = open("picture_dict.pickle", "rb")
meta = pickle.load(picture_file)
pic = meta['picture']

for ant1 in pic.keys():
    for ant2 in pic[ant1].keys():
        next_plot = np.array(pic[ant1][ant2])
        plt.plot(np.abs(next_plot[:, 0]))
        plt.show()

