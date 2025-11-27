"""Author: Andrew Martin
Creation date: 27/11/25

Script containing the functionality to produce panels for the confusion matrices, empirical cdfs and copula for the tested parametrisations.
"""

import sys
sys.path.insert(1, "../")
from common.common import get_figure_panel
from common.colormaps import (
    CMAP_probability, NORM_probability,
    CMAP_copula, NORM_copula,
    COLOR_ATL09, COLOR_Cloudnet,
    MAPPABLE_copula, MAPPABLE_probability
)
from common.handle_MI_datasets import (
    K
)
from common.handle_vcfs import vcfs_per_parametrisation, generate_confusion_matrix, generate_masks
from common.copula import BivariateCopula, patheffects

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ecdf


def plot_confusion_matrix(confusion_matrix: np.ndarray) -> (plt.Figure, plt.Axes):
    fig, ax = get_figure_panel()

    normalised_confusion_matrix = confusion_matrix / np.sum(confusion_matrix)

    ax.imshow(
        normalised_confusion_matrix,
        origin="lower",
        cmap=CMAP_probability,
        norm=NORM_probability,
    )
    for x in range(confusion_matrix.shape[1]):
        for y in range(confusion_matrix.shape[0]):
            ax.annotate(
                f"{confusion_matrix[x][y]:,}\n({normalised_confusion_matrix[x][y]:.3f})",
                xy=(x,y),
                horizontalalignment="center",
                verticalalignment="center",
                c="w",
            )

    confusion_ticks = (
        [0,1,2],
        ["nc", "pc", "tc"],
    )

    ax.set_xticks(*confusion_ticks, rotation=0)
    ax.set_yticks(*confusion_ticks, rotation=90, va="center")
    ax.set_xlabel("ATL09")
    ax.set_ylabel("Cloudnet")
    ax.set_box_aspect(1)

    return fig, ax


def plot_ecdfs(data_atl09: np.ndarray, data_cloudnet: np.ndarray, mask_non_degen: np.ndarray) -> (plt.Figure, plt.Axes):
    """Given ATL09 and Cloudnet data and a mask (all of the same shape), plot the emperical cumulative distribution functions"""
    ecdf_atl09 = ecdf(
        data_atl09[mask_non_degen].flatten()
    ).cdf
    ecdf_cloudnet = ecdf(
        data_cloudnet[mask_non_degen].flatten()
    ).cdf

    fig, ax = get_figure_panel()

    for cdf, c, label, ls in zip(
        (ecdf_atl09, ecdf_cloudnet),
        (COLOR_ATL09, COLOR_Cloudnet),
        ("ATL09", "Cloudnet"),
        ("-", (0,(5,1)))
    ):
        delta = np.ptp(cdf.quantiles) *0.05
        q = cdf.quantiles
        q = [q[0]-delta] + list(q) + [q[-1]+delta]
        u = cdf.evaluate(q)
        ax.step(
            u, q, where="post", c=c, label=label, ls=ls,
        )

    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xlabel("$u$, $v$")
    ax.set_ylabel("VCF")
    ax.set_box_aspect(1)
    ax.legend()

    return fig, ax



# load the VCF data
vcfs_per_p = vcfs_per_parametrisation()

# for each parametrisation and the associated vcfs:
# plot the confusion matrix for the vcf data
# plot the empirical cdf for the ATL09 and Cloudnet data in the (pc,pc) mask
# plot the bivariate copula density between the (pc,pc) data
for plabel, vcfs in vcfs_per_p.items():
    print(plabel)
    fig, ax = plot_confusion_matrix(
        confusion_matrix=generate_confusion_matrix(vcfs)
    )
    plt.savefig(f"{plabel}_k{K}_confusion_matrix.svg", format="svg", transparent=True)
    print("confusion matrix plot generated")

    # generate the masks that define the (pc,pc) subset of interest 
    (_, pc_atl09, _) = generate_masks(
        data = (data_atl09:=vcfs.vcf_atl09.data)
    )
    (_, pc_cloudnet, _) = generate_masks(
        data = (data_cloudnet:=vcfs.vcf_cloudnet.data)
    )
    mask_non_degenerate = pc_atl09 & pc_cloudnet
    fig, ax = plot_ecdfs(
        data_atl09=data_atl09, data_cloudnet=data_cloudnet, mask_non_degen=mask_non_degenerate
    )
    plt.savefig(f"{plabel}_k{K}_ecdf.svg", format="svg", transparent=True)
    print("ecdf plot generated")

    fig, ax = BivariateCopula.generate(
        data_X = data_atl09[mask_non_degenerate].flatten(),
        data_Y = data_cloudnet[mask_non_degenerate].flatten(),
    ).plot_density()
    plt.savefig(f"{plabel}_k{K}_copula_density.svg", format="svg", transparent=True)
    print("copula density plot generated")

# generate the colorbars for the copula density and probability
fig, cax = get_figure_panel()
plt.colorbar(mappable=MAPPABLE_copula, cax=cax)
hline = cax.axhline(1, c="k")
hline.set(path_effects = [patheffects.withTickedStroke(length=0.5, spacing=5, angle=90)])
cax.set_ylabel("copula density")
cax.set_yticks([0.5,1,1.5])
cax.yaxis.set_ticks_position('right')
cax.set_box_aspect(9)
plt.savefig(f"colorbar_copula_density.svg", format="svg", transparent=True, bbox_inches="tight")


fig, cax = get_figure_panel()
plt.colorbar(mappable=MAPPABLE_probability, cax=cax)
cax.set_ylabel("probability")
cax.set_yticks([0,1])
cax.yaxis.set_ticks_position('right')
cax.set_box_aspect(9)
plt.savefig(f"colorbar_probability.svg", format="svg", transparent=True, bbox_inches="tight")
