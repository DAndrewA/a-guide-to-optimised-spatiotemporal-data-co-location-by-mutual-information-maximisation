"""Common functions to handle loading of VCFs for given parametrisations, including p_opt at each site"""

from .common import DIR_ROOT
from .handle_MI_datasets import get_MI_with_subsetting
from .handle_sites import SITES, SITE_argument_names

import os
import xarray as xr
import numpy as np
from dataclasses import dataclass, asdict
from typing import Self
from itertools import product


dir_VCFS = os.path.join(DIR_ROOT, "vcfs_per_event")


@dataclass
class Parametrisation:
    R_km: float
    tau_s: int

    def __str__(self) -> str:
        return f"{self.R_km:3.3f}km_{self.tau_s:06}s"


@dataclass
class SiteParametrisation:
    site: str
    params: Parametrisation

    def __str__(self) -> str:
        return f"{self.site}_{self.params}"


@dataclass
class FpathsForAnalysis:
    ny_alesund: SiteParametrisation
    hyytiala: SiteParametrisation
    juelich: SiteParametrisation
    munich: SiteParametrisation

    @classmethod
    def from_parametrisation(cls, params: Parametrisation) -> Self:
        args = {
            site_arg: SiteParametrisation(site=site, params=params)
            for site, site_arg in SITE_argument_names.items()
        }
        return cls(**args)
            
    @classmethod
    def from_das_with_coords(cls, das: dict[str, xr.DataArray]) -> Self:
        fnames = {
            site: [
                print(site, R_km, tau_s),
                print(fname:=f"vcfs-per-event_{site}_{R_km:3.3f}km_{tau_s:06}s.nc"),
                fname
            ][-1]
            for site, da in das.items()
            for (R_km, tau_s) in ((float(da.R_km), int(da.tau_s)),)
        }
        return cls(**{
            site.replace("-","_"): os.path.join(dir_VCFS, fname)
            for site, fname in fnames.items()
        })

    @property
    def fpaths(self) -> dict[str, str]:
        return {
            site: os.path.join(
                dir_VCFS,   
                f"vcfs-per-event_{self.__getattribute__(site)}.nc"
            )
            for site in SITE_argument_names.values()
        }
    
    def to_DataArray(self) -> xr.DataArray:
        return xr.concat(
            [ (print(f"loading {str(fpath)}"), xr.load_dataset(str(fpath)))[-1]
              for fpath in self.fpaths.values()
            ],
            dim = "collocation_event"
        )


def site_optimised_parametrisations() -> dict[str, SiteParametrisation]:
    opt_params = dict()
    MI_full = get_MI_with_subsetting()
    for site, site_arg in SITE_argument_names.items():
        ds = MI_full.sel(site=site)
        ds = ds.isel(ds.MI.argmax(...))
        opt_params[site_arg] = SiteParametrisation(
            site=site,
            params=Parametrisation(
                R_km=float(ds.R_km),
                tau_s=int(ds.tau_s),
            )
        )
    return opt_params

_lR = [50., 500.]
_lTau = [1800, 172800]

_PARAMETRISATIONS_extremal = {
    f"P_{i}{j}": Parametrisation(R_km=R, tau_s=tau)
    for (i, R), (j, tau) in product(enumerate(_lR), enumerate(_lTau))
}

_PARAMETRISATION_literature = Parametrisation(
    R_km = 100.,
    tau_s = 10800
)

# defined as a lambda function so that it is not necesarily ran upon importing
FPATHS_by_parametrisation = lambda : (
    {
        plabel: FpathsForAnalysis.from_parametrisation(parametrisation)
        for plabel, parametrisation in _PARAMETRISATIONS_extremal.items()
    }
    | dict(P_literature = FpathsForAnalysis.from_parametrisation(_PARAMETRISATION_literature))
    | dict(P_optimal = FpathsForAnalysis(**site_optimised_parametrisations()))
)

# also defined as a lambda function so is not called upon import
vcfs_per_parametrisation = lambda: {
    plabel: fpaths.to_DataArray()
    for plabel, fpaths in FPATHS_by_parametrisation().items()
}



def generate_masks(data: np.ndarray, lower_threshold: float = 0, upper_threshold: float = 1) -> [np.ndarray, np.ndarray, np.ndarray]:
    mask_above_lower = data > lower_threshold
    mask_below_upper = data < upper_threshold
    m_lower = ~mask_above_lower
    m_middle = mask_above_lower & mask_below_upper
    m_upper = ~mask_below_upper
    return m_lower, m_middle, m_upper



def generate_confusion_matrix(vcfs: xr.Dataset) -> np.ndarray:
    """Accepts VCF datasets and computes a confusion matrix between the ATL09 and Cloudnet VCF distributions for nc, pc and tc cases.
    """
    atl09_masks = generate_masks(
        data=vcfs.vcf_atl09.data
    )
    cloudnet_masks = generate_masks(
        data=vcfs.vcf_cloudnet.data
    )
    
    confusion_matrix = np.array([
        [
            (mask_cloudnet & mask_atl09).sum()
            for mask_atl09 in atl09_masks
        ]
        for mask_cloudnet in cloudnet_masks
    ])

    return confusion_matrix

