"""Author: Andrew Martin
Creation date: 25/11/25

Script to generate the text output for the latex table describing the Cloudnet sites and the results at these locations
"""

import sys
sys.path.insert(1, "../")
from common.rho_orbits import normalised_orbital_density_degrees
from common.MI_plots import (
    DIR_MI, K, R_slice
)
from common.handle_sites import (
    SITES,
    SITE_locations,
    SITE_print_names,
)

import xarray as xr
import os
from dataclasses import dataclass
from typing import ClassVar



@dataclass
class TableRow:
    site: str
    latitude: float
    longitude: float
    orbital_density: float
    R_hat: float
    tau_hat: int
    N_events: int
    N_profiles: int
    MI: float
    MI_std: float

    def __str__(self):
        content = " & ".join((
            self.site,
            f"${self.latitude:.1f}$",
            f"${self.longitude:.1f}$",
            f"${self.orbital_density:.2f}$",
            f"${self.R_hat:.1f}$",
            f"${self.tau_hat}$",
            f"${self.N_events}$",
            "$e$".join(f"${self.N_profiles: .2e}$".split("e")),
            f"${self.MI:.3f}$ $\\pm{self.MI_std:.3f}$"
        ))
        return content + r" \\"



@dataclass
class Table:
    label: str
    caption: str
    rows: list[TableRow]
    headings: ClassVar = (
        "site",
        r"latitude (\unit{\degree N})",
        r"longitude (\unit{\degree E})",
        r"$\rho_{\text{orbits}}$",
        r"$\hat{R}$ (\unit{km})",
        r"$\hat{\tau}$ (\unit{s})",
        r"$N_{\text{events}} (\hat{\vec{p}})$",
        r"$N_{\text{profiles}} (\hat{\vec{p}})$",
        r"$\hat{\text{I}}_\text{KSG} (\hat{\vec{p}})$ (\unit{bits})",
    )

    def __str__(self):
        lines = [
            r"%t",
            r"\begin{table*}[t]",
            r"\caption{",
            f"{self.caption}",
            r"}",
            r"\label{" + f"{self.label}" + "}",
            r"\begin{tabular}{lcccccccr}",
            r"\tophline",
            " & ".join(self.headings) + r" \\",
            r"\middlehline",
            *[str(row) for row in self.rows],
            r"\bottomhline",
            r"\end{tabular}",
            r"\belowtable{} % table footnotes",
            r"\end{table*}",
        ]
        return "\n".join(lines)


caption = "\n\t".join([
    "The locations of the Cloudnet sites used in the analysis, and important results of the mutual information calculation between the ATL09 and Cloudnet VCF profiles at each site."
    r"$\rho_\text{orbits}$ represents the normalised across-track density of ICESat-2 orbits at the latitude of the Cloudnet site."
    r"$\hat{\vec{p}} = (\hat{R}, \hat{\tau})$ represents the optimised parametrisation at which the maximum mutual information, $\hat{\text{I}}_\text{KSG}(\hat{\vec{p}})$, is found."
    r"$N_\text{events}(\vec{p})$ is the number of co-location events from which data is included with a parametrisation $\vec{p}$."
    r"$N_\text{profiles}(\vec{p})$ is the number of pairwise profile comparisons made between ATL09 and Cloudnet VCF profiles across all co-location events for a given parametrisation $\vec{p}$."
])



# load in the full dataset
ds_full = xr.load_dataset(
    os.path.join(
        DIR_MI,
        "MI_merged.nc"
    )
).drop_dims("height").sel(
    dict(
        R_km=R_slice,
        K=K,
    )
)

table_rows = list()
for site in SITES:
    ds = ds_full.sel(site=site)
    maximum = ds.MI.argmax(...)
    dsmax = ds.isel(maximum)

    table_rows.append(
        tr:=TableRow(
            site = SITE_print_names[site],
            latitude=(lat:=SITE_locations[site]["lat"]),
            longitude=SITE_locations[site]["lon"],
            orbital_density=normalised_orbital_density_degrees(lat, 92.00),
            R_hat=float(dsmax.R_km),
            tau_hat=int(dsmax.tau_s),
            N_events = int(dsmax.N_events),
            N_profiles = int(dsmax.N_profiles),
            MI=float(dsmax.MI),
            MI_std=float(dsmax["std"])
        )
    )
    
    print(site, tr.R_hat * tr.orbital_density / tr.N_events)

table = Table(rows=table_rows, label="table::sites-MI", caption=caption)

with open("table-1.txt", "w") as f:
    f.write(str(table))
