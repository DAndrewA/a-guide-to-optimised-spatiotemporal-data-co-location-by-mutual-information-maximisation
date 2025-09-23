"""Author: Andrew Martin
Creation date: 22/9/25

Class definition for the areal cloud fraction homogenised data used in the analysis.
"""


from sat_val_framework.implement import HomogenisedData

import xarray as xr
import numpy as np



class CommonACFThresholds:
    """Class holding threshold values for calculating ACF values, requiring certain column cloud fractions to be met.
    """
    def __init__(self, values: np.ndarray):
        assert ( np.asarray(values) >= 0 & np.asarray(values) <= 1).all(), f"ACF threshold values must be between 0 and 1."
        self.values = np.sort(values)



class ACF(HomogenisedData):
    thresholds = None

    def assert_on_creation(self):
        assert isinstance(self.thresholds, CommonACFThresholds), f"Use an implementation of ACF with thresholds attribute set."
        assert isinstance(self.data, xr.DataArray), f"Homogenised data should be of type {xr.DataArray}, is {type(self.data)}"

        assert (
            (self.data >= 0).all()
            & (self.data <= 1).all()
        ), f"All areal cloud fractions must be within bounds [0,1]"



class ACF_10(ACF):
    thresholds = CommonACFThresholds(
        values=np.arange(0,1,0.1)
    )
