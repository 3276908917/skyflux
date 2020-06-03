import pickle
import glob
import numpy as np

"""
For each line in (argument is file path)
Delimiter is |
We assume that we only have to deal with 4 fields:
| name | right ascension in hour format | declination in arc format | flux |
"""

"""
I may have gone a little overboard with the frequency coverage;
I should definitely ask Ridhima if I should axe some of these.

GLEAMEGCAT quick guide
parameters to include in the output (downloaded to file)
name
ra
dec -40 .. -20
int_flux_76_mhz
int_flux_84_mhz
int_flux_92_mhz
int_flux_99_mhz
int_flux_107_mhz
int_flux_115_mhz
int_flux_122_mhz
int_flux_130_mhz
int_flux_143_mhz
int_flux_151_mhz
int_flux_158_mhz
int_flux_166_mhz
int_flux_174_mhz
int_flux_181_mhz
int_flux_189_mhz
int_flux_197_mhz
int_flux_204_mhz
int_flux_212_mhz
int_flux_220_mhz
int_flux_227_mhz

Alright, we are at 86633 results for dec -40 .. -20.
I may eventually simply order by 151 MHz intensity and cut at 1000 entries.
"""

# This is hard-coded to the GLEAMEGCAT format

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

        """
        We want a loop over elements from an arry. Each element describes the next frequency
        which we are extracting.
        """

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

        self.ra_angle = 15 * self.ra_hour + self.ra_minute / 4 + self.ra_second / 240

    def format_dec(self):
        remainder = self.dec
        self.dec_degree = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcminute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcsecond = float(remainder)

        self.dec_angle = self.dec_degree + self.dec_arcminute / 60 + self.dec_arcsecond / 3600

    def __str__(self):
        return "Name: " + self.name + "\nRight ascension: " + str(self.ra_angle) + \
            "\nDeclination: " + str(self.dec_angle) + \
            "\n151 MHz flux: " + self.flux_by_frq[151] + "\n"
    # we will probably want a __repr__ function so that we can see ALL fluxes associated
    # with the object.

f = open("gleam_excerpt.txt", "r")

obj_catalog = []

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

# Now print every baseline, without duplicating

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

# s(r, nu) is the Stokes I parameter.
# s = np.array([intensity goes here , 0, 0, 0])
    # and we hope that this is a 4x1 vector

"""
The following function was written by C. D. Nunhokee,
'genVisibility.py', polarizedSims, Feb 8 2019
https://github.com/Chuneeta/polarizedSims/blob/master/genVisibility.py
"""
def raddec2lm(ra0, dec0, ra, dec):# ra and dec in radians
    """
    Converts ra/dec to direction cosines l/m
    ra0  : reference/phase right ascension; type: float
    dec0 : reference/phase declination; type:float
    ra   : right ascension in radians; type:float
    dec  : declination in radians; type:float
    """
    rad2deg = lambda val: val * 180. / np.pi
    l = np.cos(dec) * np.sin(ra0 - ra)
    m = -1 * (np.sin(dec) * np.cos(dec0) - \
        np.cos(dec) * np.sin(dec0) * np.cos(ra-ra0))
    return l, m

"""
I am not so sure anymore that I really want an integral just yet.
Look at the infinitesimal: it is a solid angle.
Does that not suggest we are integrating over the whole sky?
Currently, our goal is just to pass one point source through the pipe line.
Ergo, I think it safe just to evaluate
A(r, nu) * s(r, nu) exp[-2 pi i nu b * r / c]

A(r, nu) = S^{-1} * [J(r, nu) cross J^*(r, nu)] * S

Am I evaluating this integrand for every possible base line? Would I simply sum up the integrands?
"""
