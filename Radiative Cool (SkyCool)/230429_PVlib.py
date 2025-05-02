
import numpy as np
import pandas as pd
import pvlib
from pvlib.tools import sind
from pvlib._deprecation import warn_deprecated
from pvlib.tools import _get_sample_intervals
import scipy
import scipy.constants
import warnings


skytemp= pvlib.temperature.faiman_rad(0, 20, wind_speed=1.0, ir_down=None, u0=25.0, u1=6.84, sky_view=1.0, emissivity=0.88)

print(skytemp)


