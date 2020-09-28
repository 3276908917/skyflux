import glob, os
from RIMEz import beam_models
from spin1_beam_model import cst_processing

frqs = []

data_dir = os.path.dirname(os.path.abspath(__file__)) + \
           "/skyflux/patterns/"
print("data_dir is", data_dir)

# Some files have a space after the last underscore in what I guess is
# a representation of a zero in the hundreds place. Fortunately, this
# implementation correctly ignores white space in casting the string.
file_prefix = "HERA_4.9m_E-pattern_"
file_suffix = "MHz.txt"
common_unit = 1e6 # one MHz

# generated based on magic strings above

label_start = len(data_dir + file_prefix)

pattern_list = []

for file in glob.glob(data_dir + "*"):
    print("File is:", file)
    label_end = len(file) - len(file_suffix)
    label = file[label_start:label_end]
    frq = float(label) * common_unit
    frqs.append(frq)
    print("Frq is:", frq, "\n")

    pattern_list.append(file)

model_name ="ant.h5"

processor = cst_processing.CSTDataProcessor(
    pattern_list,
    np.array(frqs),
    1, 1e-4 # magic...
)

processor.compute_spin1_harmonics()
processor.write_model_data(data_prefix, model_name)

print("Successfully output ant.h5 file. Found the following frequencies (MHz):")
print(frqs)
