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

antenna_positions = dict(pickle.load(open("ant_dict.pk", "rb")))

def baseline(ant_ID1, ant_ID2):
    return antenna_positions[ant_ID2] - antenna_positions[ant_ID1]

def list_baselines(ant_ID):
    print ("Baselines between antenna " + str(ant_ID) + " and antenna...")
    for ID in antenna_positions:
        if ant_ID != ID:
            #print(baseline(ant_ID, ID))
            print(str(ID) + ": " + str(baseline(ant_ID, ID)))

# Now print every baseline, without duplicating

# b = (u, v, w) is the
# vector representing the coordinates in meters in the plane
# of the array

# s is the Stokes I parameter?

# calculation of the celestial unit vector r
