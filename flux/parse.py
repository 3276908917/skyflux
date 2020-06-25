import pickle
import os
import numpy as np

from flux import rot
from flux import stokes

data_prefix = os.path.dirname(os.path.abspath(__file__)) + "/"
print("Searching for data files at: " + data_prefix)

# If the source catalog or antannae positions fail to load,
# we warn that we will not define the last two functions.
full_load = True

# The following section is hard-coded to the GLEAMEGCAT format,
# as downloaded by myself.
# See resources/GLEAM_guide.txt for more details.

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
        # at the same time, we convert mJy -> Jy
        for expected_frq in expected_frequencies:
            try:
                self.flux_by_frq[expected_frq] = \
                    float(line[:line.index("|")].strip()) / 1000
            except ValueError:
                print("Missing flux value for:", self.name,
                      "at frequency:", expected_frq, "MHz.")
                self.flux_by_frq[expected_frq] = np.NaN
            line = line[line.index("|") + 1:]

        try:
            self.alpha = float(line[:line.index("|")])
        except ValueError:
            print("Missing spectral index for:", self.name)
            self.alpha = np.NaN

    def format_ra(self):
        remainder = self.ra
        self.ra_hour = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_minute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_second = float(remainder)

        self.ra_angle = rot.collapse_hour(
            self.ra_hour, self.ra_minute, self.ra_second)

    def format_dec(self):
        remainder = self.dec
        self.dec_degree = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcminute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcsecond = float(remainder)

        self.dec_angle = rot.collapse_angle(
            self.dec_degree, self.dec_arcminute, self.dec_arcsecond)

    def __str__(self):
        return "Name: " + self.name + "\nRight ascension: " + str(self.ra_angle) + \
            "\nDeclination: " + str(self.dec_angle) + \
            "\n151 MHz flux: " + str(self.flux_by_frq[151]) + "\n"
    
    # we will probably want a __repr__ function so that we can see
    # ALL fluxes associated with the object.

try:
    f = open(data_prefix + "gleam_with_alpha.txt", "r")
    obj_catalog = []
    # For each line in f, the delimiter is |
    for line in f:
        obj_catalog.append(GLEAM_entry(line[1:]))
    f.close()
except FileNotFoundError:
    print("Failure to load gleam catalog.")
    full_load = False

"""
We cannot put ra0 = get_lst() in the function header. Why?
Because Python evaluates all function headers once upon first opening the script.
Consequently, the default argument will be constant: subsequent calls of this
function will use the value of get_lst() calculated when the script was opened.
"""
