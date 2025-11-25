"""Author: Andrew Martin
Creation date: 25/11/25

Script implementing the code from figure_4_5_MIs.ipynb, to generate panels for the figures showing the case study of Julich for K=10
"""

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_N_events, CMAP_N_profiles, CMAP_MI
from common.MI_plots import (
    DIR_MI,
    TEX_MI, TEX_R, TEX_tau, TEX_N_events, TEX_N_profiles,
    PLOT_ARGS, 
    K, 
    R_ticks, TAU_ticks, TAU_ticks_minor,
    R_slice,
    within_max_MI_pvalue,
    plot_data_with_maxMI_and_significance
)

import xarray as xr
import numpy as np
import os
from scipy import stats
from itertools import product
from tqdm import tqdm

import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm, Normalize

SITE = "juelich"
FIGSIZE = (12,3)



def plot_results(ds) -> (plt.Figure, list[plt.Axes]):
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True, figsize=FIGSIZE, layout="constrained")
    (ax_events, ax_profiles, ax_mi) = axs

    plot_data_with_maxMI_and_significance(
        da = ds.N_events,
        ds = ds,
        TEX_da = TEX_N_events,
        CMAP = CMAP_N_events,
        NORM_func = lambda : LogNorm(),
        ax = ax_events,
    )

    plot_data_with_maxMI_and_significance(
        da = ds.N_profiles,
        ds = ds,
        TEX_da = TEX_N_profiles,
        CMAP = CMAP_N_profiles,
        NORM_func = lambda : LogNorm(),
        ax = ax_profiles,
    )

    plot_data_with_maxMI_and_significance(
        da = ds.MI,
        ds = ds,
        TEX_da = TEX_MI,
        CMAP = CMAP_MI,
        NORM_func = lambda : Normalize(vmin=ds.MI.min(), vmax=ds.MI.max()),
        ax = ax_mi,
    )

    for ax in axs.flatten():
        ax.set_title(None)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(TEX_R)
        ax.set_ylabel(TEX_tau)
        ax.set_xticks(*R_ticks)
        ax.set_yticks(*TAU_ticks)
        ax.set_yticks(*TAU_ticks_minor, minor=True)
        ax.set_box_aspect(1)

    return fig, axs



print("Loading dataset: ", end="")
ds = xr.load_dataset(
    os.path.join(
        DIR_MI,
        "MI_merged.nc"
    )
).drop_dims("height").sel(
    dict(
        R_km=R_slice,
        site=SITE,
        K=K,
    )
)
print("success")
print(ds)
print("Computing pvalues")
ds["pvalue_MImax"] = within_max_MI_pvalue(ds, "std")
print("success")


print("Generating plot")
fig, axs = plot_results(ds)
plt.savefig(fname_out:=f"{SITE}_k{K}_results.svg", format="svg", transparent=True, bbox_inches="tight")
print(f"success, saved to {fname_out}")
