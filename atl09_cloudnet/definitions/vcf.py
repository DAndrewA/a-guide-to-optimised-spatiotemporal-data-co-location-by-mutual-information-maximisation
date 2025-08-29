"""Author: Andrew Martin
Creation date: 27/08/2025

Class definition for the Vertical cloud fraction homogenised data that is used in the analysis.
"""


from sat_val_framework.implement import HomogenisedData 

import xarray as xr
import numpy as np



class CommonGrid:
    """CommonGrid handles height values and histogram bin edges for interpolating RawData from different height grids onto a common vertical grid"""
    def __init__(self, z_values: np.ndarray, z_min:float, z_max:float):
        assert z_min < np.min(z_values), f"{z_min=:.3f} must be less than minimum supplied z value ({np.min(z_values):.3f})"
        assert z_max > np.max(z_values), f"{z_max=:.3f} must be larger than maximum supplies z value ({np.max(z_values):.3f})"
        
        self.z_values = z_values
        self.z_bounds = np.array([
            z_min, *( (z_values[:-1]+z_values[1:]) / 2 ), z_max
        ])

        self.lin_interp_z = xr.DataArray(
            data=z_values,
            dims=("height"),
            name="height"
        )
        


class VCF(HomogenisedData):
    CG: CommonGrid

    def assert_on_creation(self):
        assert isinstance(self.data, xr.DataArray), f"homogenised data should be {xr.DataArray}, is type(self.data)"
        assert (self.data.height == CG.height).all(), f"Homogenised VCF height coordinates should match common grid"
        assert (
            (self.data >= 0).all() &
            (self.data <= 1).all()
        ), f"All vertical cloud fractions must fall within bounds [0,1]"




class VCF_240m(VCF):
    CG = CommonGrid(
        z_values = np.arange(0, 12_000, 240),
        z_min = -100,
        z_max = 12_100
    )
