def generate_model:
    import os
    from RIMEz import beam_models
    from spin1_beam_model import cst_processing

    data_prefix = os.path.dirname(os.path.abspath(__file__)) + "/"
    beam_origin = data_prefix + "HERA_4.9m_E-pattern_151MHz.txt"
    model_name ="ant.h5"

    processor = cst_processing.CSTDataProcessor(
        [beam_origin, beam_origin, beam_origin],
        # hard coding to demonstrate functionality elsewhere.
        # We will fix it eventually.
        np.array([150e6, 151e6, 152e6]),
        1, 1e-4
    )

    processor.compute_spin1_harmonics()
    processor.write_model_data(data_prefix, model_name)
