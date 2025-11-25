"""Common functions to handle loading of VCFs for given parametrisations, including p_opt at each site"""
import handle_MI_datasets
from handle_sites import SITES, SITE_argument_names

import os
import xarray as xr
from dataclasses import dataclass, asdict
from typing import Self




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
        print(self.__getattribute__("ny_alesund"))
        print(self)
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
        