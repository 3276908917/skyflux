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


