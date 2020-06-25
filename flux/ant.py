"""
Utilities for calculations based on antenna positions,
such as baseline and phase factor.
"""

c = 299792458 # m / s

try:
    ant_pos = dict(pickle.load(open(data_prefix + "ant_dict.pk", "rb")))

    def baseline(ant_ID1, ant_ID2):
        """
        Calculate the baseline between antennae # @ant_ID1 and @ant_ID2
        by a simple difference of their coordinates.
        """
        return ant_pos[ant_ID2] - ant_pos[ant_ID1]

    def phase_factor(ant1, ant2, r, nu=151e6):
        """
        Calculate the phase factor in the direction @r (l, m)
            (we assume that n is of insignificant magnitude)
        and at the frequency @nu
        between two antennae whose ID #s are @ant1 and @ant2.
        When we calculate the baseline (u, v, w), we
            assume that w is of insignificant magnitude.
        """
        b = baseline(ant1, ant2)[0:2] # kill w
        br = np.dot(b, r)
        return np.exp(-2j * np.pi * nu * br / c)
    
except FileNotFoundError:
    print("Failure to load antennae data.")
