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
        line = line[line.index("|") + 1:]
        self.dec = line[:line.index("|")]
        line = line[line.index("|") + 1:]
        self.flux = line[:line.index("|")].strip()

    def __str__(self):
        return "Name: " + self.name + "\nRight ascension: " + self.ra + \
            "\nDeclination: " + self.dec + "\nWide-field flux: " + self.flux + "\n"

    #needs formatting
    def get_ra():
        return self.ra

    #needs formatting
    def get_dec():
        return self.dec

f = open("gleam_excerpt.txt", "r")

for line in f:
    print(GLEAM_entry(line[1:]))
f.close()

# Now we have a section on reading the antenna positions
    # calculate baseline vectors and all that
# is the format u, v, m?? That doesn't make any sense!
antenna_positions = dict(pickle.load(open("ant_dict.pk", "rb")))

# b = (u, v, w) is the
# vector representing the coordinates in meters in the plane
# of the array

# s is the Stokes I parameter?

# things I am currently planning to leave open
    # calculation of the celestial unit vector?
