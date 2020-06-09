import pickle
import glob
import numpy as np
import time
import astropy
import astropy.time

import rotations

"""
In principle, it should not matter which flux-by-frequency
you order by. Like, you can choose 151e6 Hz

You could also sort by int_flux (the total thing)

1 Jansky = 1000 mJy threshold would be good

get brightest source,
find distribution of source brightnesses

PARSE SPECTRAL INDICES FOR EACH OBJECT
raise a warning if an object could not find a valid spectral index!
Because, if it lacks one, we need to fit a power law.
"""

# The following section is hard-coded to the GLEAMEGCAT format

# all numbers represent MHz quantities
expected_frequencies = [76, 84, 92, 99, 107, 115, 122, 130,
                    143, 151, 158, 166, 174, 181, 189,
                    197, 204, 212, 220, 227]

class GLEAM_entry:
    def __init__(self, line):
        # Might want to redo this line later to exclude universal "GLEAM " prefix
        self.name = line[:line.index("|")]
        line = line[line.index("|") + 1:]
        
        self.ra = line[:line.index("|")]
        self.format_ra()
        line = line[line.index("|") + 1:]

        self.dec = line[:line.index("|")]
        self.format_dec()
        line = line[line.index("|") + 1:]

        self.flux_by_frq = {}

        # we extract and record fluxes according to expected_frequencies
        for expected_frq in expected_frequencies:
            self.flux_by_frq[expected_frq] = line[:line.index("|")].strip()
            line = line[line.index("|") + 1:]
            

    def format_ra(self):
        remainder = self.ra
        self.ra_hour = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_minute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_second = float(remainder)

        self.ra_angle = rotations.collapse_hour(
            self.ra_hour, self.ra_minute, self.ra_second)

    def format_dec(self):
        remainder = self.dec
        self.dec_degree = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcminute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcsecond = float(remainder)

        self.dec_angle = rotations.collapse_angle(
            self.dec_degree, self.dec_arcminute, self.dec_arcsecond)

    def __str__(self):
        return "Name: " + self.name + "\nRight ascension: " + str(self.ra_angle) + \
            "\nDeclination: " + str(self.dec_angle) + \
            "\n151 MHz flux: " + self.flux_by_frq[151] + "\n"
    # we will probably want a __repr__ function so that we can see
    # ALL fluxes associated with the object.

f = open("gleam_excerpt.txt", "r")

obj_catalog = []

# For each line in f, the delimiter is |
for line in f:
    obj_catalog.append(GLEAM_entry(line[1:]))
f.close()

# Antenna section

ant_pos = dict(pickle.load(open("ant_dict.pk", "rb")))

def baseline(ant_ID1, ant_ID2):
    """
    Calculate the baseline between antennae # @ant_ID1 and @ant_ID2
    by a simple difference of their coordinates.
    """
    return ant_pos[ant_ID2] - ant_pos[ant_ID1]

active_ants = list(ant_pos)
active_ants.sort()

S = .5 * np.array([[1, 1, 0, 0,],
                  [0, 0, 1, 1j],
                  [0, 0, 1, -1j],
                  [1, -1, 0, 0]])

c = 299792458 # m / s

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

def visibility_integrand(J, source, nu=151e6):
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([complex(I), 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    r = rotations.raddec2lm(ra, dec)
    
    J_outer = np.kron(J, np.conj(J))
    A = np.dot(np.dot(np.linalg.inv(S), J_outer), S)

    # Do I want to sum for all possible baselines?
    # No. I should make a function to evaluate the integrand
        # for a single baseline.
    phi = phase_sum(r, nu)
    return np.dot(np.dot(A, s), phi)

"""
Sum your sources, not your baselines.
Summation over different values r-hat

Integration is a summation over sources.

(I think the output should be a scalar, not a 4x1)

The graph of times and positions

A will change according to time. Equation 3 assumes that we have
just one integration time.

The change in A over time represents a drift scan,
like the Gaussians that you recently visualized.

Time you can do, you have all the information

use constant flux independent of frequency
(this is the same as saying that alpha is equal to one)
then varying frequency from 100-200 MHz
'assume it is the same beam between 100 and 200 MHz'
    beam variation is what we will eventually be handling anyway

15 m or 30 m baseline East-West

I need to plot the A matrix, to see the leakage terms.
Plot it in healpix.

(Ask for her power spectrum notebook plot for a single baseline,
if you complete all of your work.)
"""

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
