#! /usr/bin/env python

"power plot sketch"

import skyflux as sf

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

import pickle

picture_file = open("picture_dict.pickle", "rb")

meta = pickle.load(picture_file)

frq = meta['frequencies']

nu_idxs = range(len(frq))

ap = meta['ant_pos']
pic = meta['picture']

import copy
avg_pic = pic.copy() # terrible on RAM, but automatically provides a skeleton

wedge_data = []

for ant1 in pic.keys():
    for ant2 in pic[ant1].keys():
        next_plot = np.array(pic[ant1][ant2])
        plt.plot(np.abs(next_plot[:, 0]))
        plt.show()

"""
Begin scratchpad

Effective baseline is 200 MHz? Since I am running from 50-250 MHz.

We are trying to build a wedge plot. For each power, we
calculate the k proxies...

Game plan: let us have some universal array
	wedge_data
and dump triples (k_orth, k_parr, p_p) into it,
	then we can use matplotlib to get a visual

"""


