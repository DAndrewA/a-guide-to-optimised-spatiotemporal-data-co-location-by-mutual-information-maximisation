"""Author: Andrew Martin
Creation date: 25/11/25

Script implementing the code from figure_4_5_MIs.ipynb, to generate panels for the figures showing the mutual information surfaces at all sites for K=10
"""

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_N_events, CMAP_N_profiles, CMAP_MI
from common.handle_sites import SITES
from common.handle_MI_datasets import get_MI_with_subsetting, K
from common.MI_plots import (
    TEX_MI, TEX_R, TEX_tau,
    PLOT_ARGS, FIGSIZE,
    R_ticks, TAU_ticks, TAU_ticks_minor,
    within_max_MI_pvalue,
    plot_data_with_maxMI_and_significance,
)

import xarray as xr
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerLine2D

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
ds_full = get_MI_with_subsetting()
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



# legend
def legend_handles():
    """Return legend handles for a vertical, horizontal dashed line, and hatching"""

    def update_prop(handle, orig):
        handle.update_from(orig)
        x,y = handle.get_data()
        handle.set_data([np.mean(x)]*2, [0, 2*y[0]])

    handles = [
        (vertical_line:=Line2D([0.5,0.5],[0,1], c="k", ls="--", lw=1, label=r"$\hat{R}$")),
        Line2D([],[], c="k", ls="--", lw=1, label=r"$\hat{\tau}$"),
        Patch(facecolor="none", hatch="xxx", label=r"possible $\hat{\boldsymbol{p}}$ candidate") 
    ]
    handler_map = {
        vertical_line: HandlerLine2D(update_func=update_prop)
    }
    return handles, handler_map


fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")
handles, handler_map = legend_handles()
legend = plt.legend(handles=handles, handler_map=handler_map, loc="upper left", handlelength=1, handleheight=1)
ax.axis("off")
plt.savefig("legend.svg", format="svg", transparent=True, bbox_inches="tight")
