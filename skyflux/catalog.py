import numpy as np

from skyflux import rot

import warnings as w
import os

# The following section is hard-coded to the GLEAMEGCAT format,
# as downloaded by myself.
# See resources/GLEAM_guide.txt for more details.

# all numbers represent MHz quantities
expected_frequencies = [76, 84, 92, 99, 107, 115, 122, 130,
                    143, 151, 158, 166, 174, 181, 189,
                    197, 204, 212, 220, 227]

class GLEAM_entry:
    def __init__(self, line):
        """
        Initialize GLEAM object with a line from the catalog
        excerpt downloaded as a text file. We expect all of
        the following, in this order, on a separate line
        associated with a single GLEAM object:
            name
            right ascension (hour, minute, second)
            declination (degree, arcminute, arcsecond)
            integrated flux value for each of the frequencies in
                expected_frequencies
            spectral index, alpha
        """
        # Might want to redo this line later to
        # exclude universal "GLEAM " prefix
        self.name = line[:line.index("|")]
        line = line[line.index("|") + 1:]
        
        self.ra = line[:line.index("|")]
        self._format_ra()
        line = line[line.index("|") + 1:]

        self.dec = line[:line.index("|")]
        self._format_dec()
        line = line[line.index("|") + 1:]

        self.flux_by_frq = {}

        # we extract and record fluxes according to expected_frequencies
        # at the same time, we convert mJy -> Jy
        for expected_frq in expected_frequencies:
            try:
                self.flux_by_frq[expected_frq] = \
                    float(line[:line.index("|")].strip()) / 1000
            except ValueError:
                warning = "Missing flux value for: " + self.name + \
                      " at frequency: " + str(expected_frq) + " MHz."
                w.warn(warning)
                self.flux_by_frq[expected_frq] = np.NaN
            line = line[line.index("|") + 1:]

        try:
            self.alpha = float(line[:line.index("|")])
        except ValueError:
            warning = "Missing spectral index for: " + self.name
            w.warn(warning)
            self.alpha = np.NaN

    def _format_ra(self):
        """
        self.ra is a string which describes the right-ascension
            of the source in the format 'HH-mm-ss'.
        This internal function breaks that string up into
            three explicitly floating-point instance variables:
                ra_hour, ra_minute, and ra_second
        Finally, it allocates the instance variable ra_angle
            to represent a full conversion of the original
            right-ascension into a single, centralized, floating-point
            angle in degrees.
        """
        remainder = self.ra
        self.ra_hour = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_minute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.ra_second = float(remainder)

        self.ra_angle = rot.collapse_hour(
            self.ra_hour, self.ra_minute, self.ra_second)

    def _format_dec(self):
        """
        self.dec is a string which describes the declination
            of the source in the format
            '(arc)degree-arcminute-arcsecond'.
        This internal function breaks that string up into
            three explicitly floating-point instance variables:
                dec_degree, dec_arcminute, and dec_arcsecond
        Finally, it allocates the instance variable dec_angle
            to represent a full conversion of the original
            declination into a single, centralized, floating-point
            angle in degrees.
        """
        remainder = self.dec
        self.dec_degree = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcminute = float(remainder[:remainder.index(" ")])

        remainder = remainder[remainder.index(" ") + 1:]
        self.dec_arcsecond = float(remainder)

        self.dec_angle = rot.collapse_angle(
            self.dec_degree, self.dec_arcminute, self.dec_arcsecond)

    def __str__(self):
        return "Name: " + self.name + \
            "\nRight ascension: " + str(self.ra_angle) + \
            "\nDeclination: " + str(self.dec_angle) + \
            "\n151 MHz flux: " + str(self.flux_by_frq[151]) + "\n"
    
    # we will probably want a __repr__ function so that we can see
    # ALL fluxes associated with the object.

data_prefix = os.path.dirname(os.path.abspath(__file__)) + "/"

try:
    srcs = np.load(data_prefix + 'catalog.npy', allow_pickle=True)
except FileNotFoundError:
    print("Failure to load GLEAM object array.")
    
def idxs_to_objs(indices):
    objs = []
    for i in indices:
        objs.append(srcs[i])
    return np.array(objs)

