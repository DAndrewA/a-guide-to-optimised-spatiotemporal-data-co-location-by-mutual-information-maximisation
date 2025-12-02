"""Author: Andrew Martin
Creation date: 25/11/25

Script implementing the code from figure_4_5_MIs.ipynb, to generate panels for the figures showing the case study of Julich for K=10
"""

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_N_events, CMAP_N_profiles, CMAP_MI
from common.handle_MI_datasets import (
    get_MI_with_subsetting,
    K,
)
from common.MI_plots import (
    TEX_MI, TEX_R, TEX_tau, TEX_N_events, TEX_N_profiles,
    PLOT_ARGS, FIGSIZE,
    R_ticks, TAU_ticks, TAU_ticks_minor,
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
ds = get_MI_with_subsetting(site=SITE).load()
print("success")
print(ds)

print(f"argmin: ", ds.isel(ds.MI.argmin(...)))
print(f"argmax: ", ds.isel(ds.MI.argmax(...)))

print("Computing pvalues")
ds["pvalue_MImax"] = within_max_MI_pvalue(ds, "std")
print("success")


print("Generating plots")
for da, TEX_da, CMAP, NORM_func, savename in (
    (ds.N_events, TEX_N_events, CMAP_N_events, (lambda : LogNorm()), "n_events"),
    (ds.N_profiles, TEX_N_profiles, CMAP_N_profiles, (lambda: LogNorm()), "n_profiles"),
    (ds.MI, TEX_MI, CMAP_MI, (lambda : Normalize(vmin=ds.MI.min(), vmax=ds.MI.max())), "mutual_information"),
):
    fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")
    
    plot_data_with_maxMI_and_significance(
        da = da,
        ds = ds,
        TEX_da = TEX_da,
        CMAP = CMAP,
        NORM_func = NORM_func,
        ax = ax
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

    plt.savefig(fname_out:=f"{SITE}_k{K}_{savename}.svg", format="svg", transparent=True, bbox_inches="tight")
print(f"success, saved to {fname_out}")
