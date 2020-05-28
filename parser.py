import pickle
import glob

"""
For each line in (argument is file path)
Delimiter is |
We assume that we only have to deal with 4 fields:
| name | right ascension in hour format | declination in arc format | flux |

I do not know how to open a file
We want a Python trim function, because,
    when we try to read from the flux field,
    we will be ingesting a ton of white space.
"""

# This is hard-coded to the GLEAMEGCAT format

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

        self.flux = line[:line.index("|")].strip()

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
            "\nDeclination: " + str(self.dec_angle) + "\nWide-field flux: " + self.flux + "\n"

f = open("gleam_excerpt.txt", "r")

for line in f:
    print(GLEAM_entry(line[1:]))
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

# s is the Stokes I parameter. I am not sure how to get that vector out of the integral.
# I need to calculate the celestial unit vector r, but the textbook is a little over my head.
