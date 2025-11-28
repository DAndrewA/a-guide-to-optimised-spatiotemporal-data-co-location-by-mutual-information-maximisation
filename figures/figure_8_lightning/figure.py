"""Author: Andrew Martin
Creation date: 27/11/25

Script that creates the VCF mean bias profile plots.
"""

import sys
sys.path.insert(1, "../")
from common.colormaps import (
    CMAP_probability, NORM_probability,
    MAPPABLE_probability, LogNorm,
    CMAP_plabels
)
from common.handle_vcfs import vcfs_per_parametrisation, PARAMETRISATION_print_names

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from xhistogram.xarray import histogram as xhist
from dataclasses import dataclass

PANEL_SIZE = (4,5)


BINS_bias = np.linspace(-1,1,50) 
BINS_logratio = np.linspace(-3,3,60)
NORM_bias = LogNorm(vmax=np.power(10,1.5), vmin=np.power(10,-2.5), clip=True)

YBOUNDS = [0, 11_500] 
YTICKS_packed = (
    [0, 2_000, 4_000, 6_000, 8_000, 10_000],
    [0,2,4,6,8,10]
)


# class handling how the parametrisations are plotted
@dataclass
class ParametrisationPlotArgs:
    marker: str | None
    ls: str | None
    lc: tuple


def dvcfs(vcfs: xr.Dataset, threshold:float=0.) -> xr.DataArray:
    dvcfs = xhist(
        (vcfs.vcf_atl09 - vcfs.vcf_cloudnet).rename("dvcf")
            .where(
                (vcfs.vcf_atl09 > threshold)
                & (vcfs.vcf_cloudnet > threshold)
            ),
        bins = [BINS_bias],
        dim=["collocation_event"],
        density=True,
    )
    # normalise the dvcfs
    dvcfs = dvcfs / dvcfs.integrate(coord="dvcf_bin")
    return dvcfs


def logratios(vcfs: xr.Dataset) -> xr.DataArray:
    return xhist(
        np.log10(vcfs.vcf_atl09 / vcfs.vcf_cloudnet).rename("logratio"),
        bins = [BINS_logratio],
        dim=["collocation_event"],
        density=True,
    )



vcfs_per_p = vcfs_per_parametrisation()

del vcfs_per_p["P_01"], vcfs_per_p["P_10"]

dvcfs_per_p = {
    plabel: dvcfs(vcfs)
    for plabel, vcfs in vcfs_per_p.items()
}


#df = 0.15
df = 1 # when CMAP_plabels=cm.batlowS
plot_args_by_p = {
    plabel: ParametrisationPlotArgs(
        marker = m,
        ls = ls,
        lc = lc
    )
    for plabel, m, ls, lc in zip(
        dvcfs_per_p.keys(),
        ("o", "s", "x", None), # markers
        ("-.", ":", "--", "-"), # linestyles
        (CMAP_plabels(0), CMAP_plabels(df), CMAP_plabels(2*df), CMAP_plabels(3*df))  # colors
    )
}


def expected_bias(dvcf: xr.DataArray) -> xr.DataArray:
    """Compute the expected VCF difference"""
    return (dvcf * dvcf.dvcf_bin).integrate(coord="dvcf_bin")

def variance_dvcf(dvcf: xr.DataArray) -> xr.DataArray:
    """Compute the variance of the VCF difference distribution"""
    return (
        (dvcf * np.power(dvcf.dvcf_bin,2)).integrate(coord="dvcf_bin") 
        - np.power(expected_bias(dvcf),2)
    )

means_by_p = {
    plabel: expected_bias(dvcf)
    for plabel, dvcf in dvcfs_per_p.items()
}
vars_by_p = {
    plabel: variance_dvcf(dvcf)
    for plabel, dvcf in dvcfs_per_p.items()
}

def plot_distribution(plabel: str, dvcf: xr.DataArray, expected_bias: xr.DataArray) -> (plt.Figure, plt.Axes):
    fig, ax = plt.subplots(1,1, figsize=PANEL_SIZE, layout="constrained")
    plot_args = plot_args_by_p[plabel]

    dvcf.plot(
        y="height",
        cmap=CMAP_probability,
        norm=NORM_bias,
        ax=ax,
        add_colorbar=False
    )

    expected_bias.plot(
        y="height",
        ax=ax,
        lw=2,
        ls = plot_args.ls,
        c = plot_args.lc,
        marker = plot_args.marker,
        markevery=5,
    )

    ax.axvline(0, ls=(12,(10,7)), lw=1, c="w")

    ax.set_yticks(*YTICKS_packed)
    ax.set_ylim(YBOUNDS)
    ax.set_xticks([-1,0,1])
    ax.set_xlim([-1,1])

    ax.set_ylabel(r"$z$ (km)")
    ax.set_xlabel(r"$\nu$")
    ax.set_box_aspect(2.5)

    ax.set_title(PARAMETRISATION_print_names[plabel], ha="center")

    return fig, ax


def plot_expected_differences(means_by_p: dict[str, xr.DataArray]) -> (plt.Figure, plt.Axes):
    fig, ax = plt.subplots(1,1, figsize=PANEL_SIZE, layout="constrained")
    ax.set_facecolor(CMAP_probability(0))

    ax.axvline(0, ls=(12,(10,7)), lw=1, c="w")

    for i, (plabel, expected_dvcf) in enumerate(means_by_p.items()):
        plot_args = plot_args_by_p[plabel]

        expected_dvcf.plot(
            y="height",
            ax=ax,
            lw=2,
            ls = plot_args.ls,
            c = plot_args.lc,
            marker = plot_args.marker,
            markevery=(2*i,8),
        )

    ax.set_yticks(*YTICKS_packed)
    ax.set_ylim(YBOUNDS)
    #ax.set_xticks([-1,0,1])
    #ax.set_xlim([-1,1])

    ax.set_ylabel(r"$z$ (km)")
    ax.set_xlabel(r"$\mathbb{E}[\nu\,|\,z]$")
    ax.set_box_aspect(2.5)

    return fig, ax

def plot_dvcf_variances(variances_by_p: dict[str, xr.DataArray]) -> (plt.Figure, plt.Axes):
    fig, ax = plt.subplots(1,1, figsize=PANEL_SIZE, layout="constrained")
    ax.set_facecolor(CMAP_probability(0))

    for i, (plabel, variance) in enumerate(variances_by_p.items()):
        plot_args = plot_args_by_p[plabel]

        variance.plot(
            y="height",
            ax=ax,
            lw=2,
            ls = plot_args.ls,
            c = plot_args.lc,
            marker = plot_args.marker,
            markevery=(2*i,8),
            label=PARAMETRISATION_print_names[plabel],
        )

    ax.set_yticks(*YTICKS_packed)
    ax.set_ylim(YBOUNDS)
    ax.set_xlim([0,0.25])

    ax.set_ylabel(r"$z$ (km)")
    ax.set_xlabel(r"$\text{Var}[\nu \, | \, z]$")
    ax.set_box_aspect(2.5)

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    return fig, ax


# plot the vcf difference distributions
for plabel, dvcf in dvcfs_per_p.items():
    fig, ax = plot_distribution(
        plabel = plabel,
        dvcf = dvcf,
        expected_bias = means_by_p[plabel]
    )
    plt.savefig(f"{plabel}_bias_distribution.svg", format="svg", transparent=True, bbox_inches="tight")
    
fig, ax = plot_expected_differences(means_by_p=means_by_p)
plt.savefig("expected_bias_distribution.svg", format="svg", transparent=True, bbox_inches="tight")

fig, ax = plot_dvcf_variances(variances_by_p=vars_by_p)
plt.savefig("bias_variance_distributions.svg", format="svg", transparent=True, bbox_inches="tight")
