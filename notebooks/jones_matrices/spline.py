import numpy as np
import RIMEz

spline_function = RIMEz.beam_models.model_data_to_spline_beam_func('./HERA_4.9m_E-pattern_151MHz.txt', np.array([151e6]))
# What does the argument i do?


# remember that the Jones matrix depends on \hat{r} and \nu, the frequency (of what?)

# Now I just need to find what the actual arguments are...

# Are we using radians or degrees?
print(spline_function(1, 151e6, 30, 30))

