print("\nStarting script to generate spin-1 harmonics " + \
      "for the beam files stored in ./patterns/")

import glob, os
import numpy as np
from RIMEz import beam_models
from spin1_beam_model import cst_processing

frqs = []

# split up the pathing components to facilitate maintenance

# the following line cannot be replicated in any shell
here = os.path.dirname(os.path.abspath(__file__)) + "/"
data_dir = here + "patterns"
output_model_name = "HERA_spin1_harmonics"

# Some files have a space after the last underscore in what I guess is
# a representation of a zero in the hundreds place. Fortunately, this
# implementation correctly ignores white space in casting the string.
file_prefix = "/HERA_4.9m_E-pattern_"
file_suffix = "MHz.txt"
common_unit = 1e6 # one MHz

# generated based on magic strings above

label_start = len(data_dir + file_prefix)

pattern_list = []

print("\nAll header variables processed. Walking file directory...")

for file in glob.glob(data_dir + "/*"):
    label_end = len(file) - len(file_suffix)
    label = file[label_start:label_end]
    frq = float(label) * common_unit
    frqs.append(frq)

    pattern_list.append(file)

print("\nWalk complete. Generating CST data processor...")

processor = cst_processing.CSTDataProcessor(
    pattern_list,
    np.array(frqs),
    1, 1e-4 #! magic...
)

print("\nCST data processor generated. Computing spin harmonics...")

processor.compute_spin1_harmonics()

# os.path.join sucks, but CST processor uses it. Here's a hack:
processor.write_model_data("", here + output_model_name)

print("\nSuccessfully wrote spin harmonics to disk." + \
      "\nFound the following frequencies (MHz):")
print(frqs)
