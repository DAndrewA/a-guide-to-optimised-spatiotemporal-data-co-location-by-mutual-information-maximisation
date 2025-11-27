"""Common functions to handle loading of VCFs for given parametrisations, including p_opt at each site"""

from handle_MI_datasets import get_MI_with_subsetting, K, R_slice
from handle_sites import SITES, SITE_argument_names

import os
import xarray as xr
from dataclasses import dataclass, asdict
from typing import Self
from itertools import product


SCRATCH = "/work/scratch-nopw2/eeasm"
dir_VCFS = os.path.join(SCRATCH, "vcfs_per_event")


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
    MI_full = get_MI_with_subsetting(K=K, R=R_slice)
    for site, site_arg in SITE_argument_names.items():
        ds = MI_full.sel(site)
        ds = ds.isel(ds.MI.argmax(...))
        opt_params[site] = SiteParametrisation(
            site=site,
            params=Parametrisation(
                R_km=float(ds.R_km),
                tau_s=float(ds.tau_s),
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

FPATHS_by_parametrisation = (
    {
        plabel, FpathsForAnalysis.from_parametrisation(parametrisation)
        for plabel, parametrisation in _PARAMETRISATIONS_extremal.items()
    }
    | dict(P_literature = FpathsForAnalysis.from_parametrisation(_PARAMETRISATION_literature))
    | dict(P_optimal = FpathsForAnalysis(**site_optimised_parametrisations()))
)
