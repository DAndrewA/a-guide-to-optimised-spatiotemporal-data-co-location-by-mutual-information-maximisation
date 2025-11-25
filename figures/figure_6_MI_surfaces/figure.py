"""Author: Andrew Martin
Creation date: 25/11/25

Script implementing the code from figure_4_5_MIs.ipynb, to generate panels for the figures showing the mutual information surfaces at all sites for K=10
"""

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_N_events, CMAP_N_profiles, CMAP_MI
from common.handle_sites import SITES
from common.MI_plots import (
    DIR_MI,
    TEX_MI, TEX_R, TEX_tau,
    PLOT_ARGS, FIGSIZE,
    K, 
    R_ticks, TAU_ticks, TAU_ticks_minor,
    R_slice,
    within_max_MI_pvalue,
    plot_data_with_maxMI_and_significance,
)

import xarray as xr
import numpy as np
import os

import matplotlib.pyplot as plt

from matplotlib.colors import Normalize



def plot_mutual_information_panel(ds: xr.Dataset) -> (plt.Figure, list[plt.Axes]):
    fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")
    
    plot_data_with_maxMI_and_significance(
        da = ds.MI,
        ds = ds,
        TEX_da = TEX_MI,
        CMAP = CMAP_MI,
        NORM_func = lambda : Normalize(vmin=ds.MI.min(), vmax=ds.MI.max()),
        ax = ax,
    )

    ax.set_title(None)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(TEX_R)
    ax.set_ylabel(TEX_tau)
    ax.set_xticks(*R_ticks)
    ax.set_yticks(*TAU_ticks)
    ax.set_yticks(*TAU_ticks_minor, minor=True)
    ax.set_box_aspect(1)

    return fig, ax


print("Loading dataset: ", end="")
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
print("success")
print(ds_full)
for site in SITES:
    print(f"CREATING PANEL FOR {site}")
    ds = ds_full.sel(site=site)
    print("Computing pvalues")
    ds["pvalue_MImax"] = within_max_MI_pvalue(ds, "std")
    print("success")

    print("Generating plot")
    fig, axs = plot_mutual_information_panel(ds)
    plt.savefig(fname_out:=f"{site}_k{K}_mutual_information.svg", format="svg", transparent=True, bbox_inches="tight")
    plt.clf()
    print(f"success, saved to {fname_out}")
