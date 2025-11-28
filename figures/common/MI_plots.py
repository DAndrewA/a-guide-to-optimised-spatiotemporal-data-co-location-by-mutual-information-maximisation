"""Author: Andrew Martin
Creation date: 25/11/25

Script containing common functionality between the plots for figures 4(5) and 5(6). (brackets indicate the new figure number after the inclusion of figure 3 -- synthetic data regime examples)
"""

import os
import numpy as np
import xarray as xr
from matplotlib.colors import Normalize 
from scipy import stats
from tqdm import tqdm
from itertools import product
from typing import Callable

from .colormaps import CMAP_MI
from .common import PANEL_SIZE as FIGSIZE

PLOT_ARGS = lambda : dict(
    x = "R_km",
    y = "tau_s",
)

TEX_MI = r"$\hat{\text{I}}_\text{KSG} ( \boldsymbol{p} )$ (nats)"
TEX_N_events = r"$N_\text{events}$"
TEX_N_profiles = r"$N_\text{profiles}$"
TEX_R = r"$R$ (km)"
TEX_tau = r"$\tau$ (hours)"

R_ticks = (
    _temp:=[50, 100, 500],
    _temp
)
del _temp

TAU_ticks = (
    [600, 3600, 3600*3, 3600*12, 3600*24, 3600*24*2],
    #["10 minutes", "1 hour", "3 hours", "12 hours", "1 day", "2 days"]
    [r"$\frac{1}{6}$", "1", "3", "12", "24", "48"] # units of hours
)
TAU_ticks_minor = (
    mticks:=[
        *np.arange(600, 3600, 600),
        *np.arange(3600, 3600*12, 3600),
        *np.arange(3600*12, 3600*24*2+1, 3600*12)
    ],
    [None]*len(mticks)
)

def within_max_MI_pvalue(ds: xr.Dataset, std_var: str = "std") -> xr.DataArray:
    MI_max_pvalues = xr.zeros_like(ds.MI, dtype=float)
    argmax = ds.MI.argmax(dim=["R_km","tau_s"])
    ds_argmax = ds.isel(argmax)
    for R_km, tau_s in tqdm(product(ds.R_km, ds.tau_s)):
        current_sel = dict(R_km=R_km, tau_s=tau_s) 
        da = ds.sel(current_sel)
        pvalue = stats.ttest_ind_from_stats(
            mean1 = float(ds_argmax.MI),
            std1 = float(ds_argmax[std_var]),
            nobs1 = int(ds_argmax.n_splits_std),
            mean2 = float(da.MI),
            std2 = float(da[std_var]),
            nobs2 = int(da.n_splits_std),
            equal_var=False,
            alternative='greater'
        ).pvalue
        MI_max_pvalues.loc[current_sel] = pvalue
    return MI_max_pvalues


def plot_data_with_maxMI_and_significance(
    da: xr.DataArray, 
    ds: xr.Dataset, 
    TEX_da: str,
    CMAP,
    NORM_func: Callable[[],Normalize],
    ax,
):
    # Plot the data given in the data array da in three steps:
    # 1. plot the data with its natural normalising
    # 2. Overplot hatching across the whole panel
    # 3. plot the data array again, masked for where the data is significant
    NORM_MI = Normalize(vmax = ds.MI.max(), vmin = ds.MI.min())

    da.rename(TEX_da).plot(
        **PLOT_ARGS(),
        norm=NORM_func(),
        cmap = CMAP,
        ax = ax,
    )

    ax.fill_between(
        x = [0,1],
        y1=0, y2=1,
        transform = ax.transAxes,
        hatch = "xxxxx",
        color = "none",
        edgecolor = "k",
        lw = 1,
    )

    da.rename(TEX_da).where(ds.pvalue_MImax <= 0.05).plot(
        **PLOT_ARGS(),
        norm=NORM_func(),
        cmap = CMAP,
        ax = ax,
        add_colorbar=False
    )

    # plot the maximum MI value using dashwed lines
    maximum = ds.MI.argmax(...)
    ax.axvline(
        ds.MI.isel(maximum).R_km,
        ls="--", c="k", lw=1
    )
    ax.axhline(
        ds.MI.isel(maximum).tau_s,
        ls="--", c="k", lw=1
    )
