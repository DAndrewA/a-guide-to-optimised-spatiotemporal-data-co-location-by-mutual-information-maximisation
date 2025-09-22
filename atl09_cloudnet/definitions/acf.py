"""Author: Andrew Martin
Creation date: 22/9/25

Class definition for the areal cloud fraction homogenised data used in the analysis.
"""


from sat_val_framework.implement import HomogenisedData

import xarray as xr
import numpy as np



class ACF(HomogenisedData):
    """Homogenised data for the areal cloud fraction from a given data source.
    This is defined as the fraction of vertical profiles included in the analysis that contain any amount of cloud.
    This is simply stored as a float to the data field, that must be between 0 and 1.
    """

    def assert_on_creation(self):
        assert isinstance(self.data, float), f"homogenised data should be of type float, is type {self.data}".
        assert self.data >= 0 and self.data <= 1, f"self.data should be bounded in the interval [0,1], has value {self.data}"

