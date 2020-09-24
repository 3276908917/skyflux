frqs = []

def generate_model():
    global frqs

    import glob, os
    from RIMEz import beam_models
    from spin1_beam_model import cst_processing

    data_prefix = os.path.dirname(os.path.abspath(__file__)) + \
                  "HERA_4.9m_E-pattern_ "
    data_suffix = "MHz.txt"
    common_unit = 1e6 # one MHz

    # generated based on magic strings above

    label_start = len(data_prefix)

    pattern_list = []
    
    for file in glob.glob("./patterns/*"):
        label_end = len(file) - len(data_suffix)
        label = file[label_start:label_end]
        frq = float(label) * common_unit
        frqs.append(frq)

        pattern_list.append(file)
    
    model_name ="ant.h5"

    processor = cst_processing.CSTDataProcessor(
        pattern_list,
        np.array(frqs),
        1, 1e-4 # magic...
    )

    processor.compute_spin1_harmonics()
    processor.write_model_data(data_prefix, model_name)
