import pickle
import glob
import numpy as np
import time
import astropy
import astropy.time

import rotations

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

"""
b = (u, v, w) is the
vector representing the coordinates in meters in the plane
of the array
"""

def baseline(ant_ID1, ant_ID2):
    """
    Calculate the baseline between antennae # @ant_ID1 and @ant_ID2
    by a simple difference of their coordinates.
    """
    return ant_pos[ant_ID2] - ant_pos[ant_ID1]

def list_baselines(ant_ID):
    """
    Print every baseline that features antenna # @ant_ID
    """
    print ("Baselines between antenna " + str(ant_ID) + " and antenna...")
    for ID in ant_pos:
        if ant_ID != ID:
            print(str(ID) + ": " + str(baseline(ant_ID, ID)))

active_ants = list(ant_pos)
active_ants.sort()

def all_baselines():
    """
    Print every baseline, without duplicating.

    To-do: eventually I am going to want to figure out how to propagate
        the particular order (in which I am subtracting coordinates)
        out to the integral that I am taking.
    """
    for i in range(len(active_ants)):
        ID1 = active_ants[i]
        for j in range(i + 1, len(active_ants[i + 1:])):
            ID2 = active_ants[j]
            print("Baseline between antennae " + str(ID1) + \
                  " and " + str(ID2) + " = " + str(baseline(ID1, ID2)))

"""
The following function was written by C. D. Nunhokee,
'genVisibility.py', polarizedSims, Feb 8 2019
https://github.com/Chuneeta/polarizedSims/blob/master/genVisibility.py
"""
# Unfortunately, the parameter order has been switched wrt the citation.
def raddec2lm(ra, dec, ra0=None, dec0=hera_lat): # ra and dec in radians
    """
    Converts ra/dec to direction cosines l/m
    ra0  : reference/phase right ascension; type: float
    dec0 : reference/phase declination; type:float
    ra   : right ascension in radians; type:float
    dec  : declination in radians; type:float
    """
    # See note at the end about default arguments.
    if ra0 is None:
        ra0 = get_lst()

    l = np.cos(dec) * np.sin(ra0 - ra)
    m = -1 * (np.sin(dec) * np.cos(dec0) - \
        np.cos(dec) * np.sin(dec0) * np.cos(ra-ra0))
    return l, m

"""
I am not so sure that I really want an integral just yet.
Look at the infinitesimal: it is a solid angle.
Does that not suggest we are integrating over the whole sky?
Currently, our goal is just to pass one point source through the pipe line.
"""

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

def phase_sum(r, nu=151e6):
    """
    add together the phase_factors for all
    possible baselines without duplicates
    """
    total = complex(0)

    for i in range(len(active_ants)):
        ID1 = active_ants[i]
        for j in range(i + 1, len(active_ants[i + 1:])):
            ID2 = active_ants[j]
            total += phase_factor(ID1, ID2, r, nu)
            
    return total

def visibility_integrand(J, source, nu=151e6):
    I = source.flux_by_frq[nu / 1e6]
    s = np.array([I, 0, 0, 0])

    ra = np.radians(source.ra_angle)
    dec = np.radians(source.dec_angle)
    r = raddec2lm(ra, dec)
    
    J_outer = np.kron(J, np.conj(J))
    A = np.dot(np.dot(np.linalg.inv(S), J_outer), S)

    # Do I want to sum for all possible baselines?
    # If not: should I make a function to evaluate the integrand for a single baseline?
    return np.dot(np.dot(A, s), phase_sum(r, nu))

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
