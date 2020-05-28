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

f = open("gleam_excerpt.txt", "r")
for line in f:
    line = line[1:]
    line.index("|")
    print(line.strip())
f.close()

# Now we have a section on reading the antenna positions
    # calculate baseline vectors and all that
