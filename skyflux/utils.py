"""
Miscellaneous materials that are more universally applicable than
the contents of demo.py
"""

import numpy as np
import pickle

def load_saves(filepath):
    """
    Return a dictionary containing the arrays saved to
    the .npz file at @filename

    Do not include the file ending in the argument.
    """
    a = np.load(filepath, allow_pickle=True)
    return dict(zip(("{}".format(k) for k in a), (a[k] for k in a)))
    
def pickle_dict(dict_, label):
    """
    Pickles @dict_ with a file name based on @label.
    While this is a generic routine and will pickle any dictionary,
        the intent is solely for use in conjunction with the
        package_wedge routine.
    """
    with open(label + '.pickle', 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
