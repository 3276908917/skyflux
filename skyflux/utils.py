"""
Miscellaneous materials that are more universally applicable than
the contents of demo.py
"""

import numpy as np

def load_saves(filename):
    """
    Return a dictionary containing the arrays saved to
    the .npz file at @filename

    Do not include the file ending in the argument.
    """
    a = np.load(filename + '.npz', allow_pickle=True)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))
