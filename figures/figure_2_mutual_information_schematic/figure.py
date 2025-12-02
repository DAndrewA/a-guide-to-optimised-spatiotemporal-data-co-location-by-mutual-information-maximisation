"""Author: Andrew Martin
Creation date: 2/12/25

Script to produce plot showing the independent and dependent data examples from two polynomial probability distributions.
"""

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_probability, MAPPABLE_probability
from common.synthetic_data import gen_xy_samples_mixture
from common.handle_MI_datasets import K

sys.path.insert(1, "../../")
from atl09_cloudnet.holmes_et_al_2019.python_interface import MIEstimate

a,b,c = 1, 0.3, 0.05
BINS = 100

p, FIGSIZE = (0.01, (6.8, 4))
COLOR_hist = CMAP_probability(0.7)

RHOSTAR = 0.99
N_samples = 50_000

(
    fig, 
    (
        (NOO, axmx1, axmx2, NO),
        (axmy, axh1, axh2, cax)
    )
) = plt.subplots(
    2,4, 
    width_ratios=[b,a,a,c], 
    height_ratios=[b,a], 
    sharey="row", sharex="col", 
    layout="constrained", 
    figsize=FIGSIZE,
    gridspec_kw={'hspace': p+0.1, 'wspace': p}
)

# generate samples with low mutual information (indpendent)
(_,_), (x_indep, y_indep) = gen_xy_samples_mixture(
    n_samples=N_samples,
    rhostar=RHOSTAR,
    mixture=0,
)
y_indep = 1 - y_indep

h1 = axh1.hist2d(
    x_indep, y_indep,
    bins=BINS,
    density=True,
    cmap=CMAP_probability,
)
axh1.set_yticks([])


# generate and plot dependent data
(_,_), (x_dep, y_dep) = gen_xy_samples_mixture(
    n_samples=N_samples,
    rhostar=RHOSTAR,
    mixture=1
)
y_dep = 1-y_dep

h2 = axh2.hist2d(
    x_dep, y_dep,
    bins=BINS,
    density=True,
    cmap=CMAP_probability,
)

# colorbar
plt.colorbar(MAPPABLE_probability, cax=cax)
cax.set_ylabel("Normalised probability density")
cax.set_yticks([0,1])
cax.yaxis.set_ticks_position("left")


# Y marginal distribution: only a singular plot
axmy.invert_xaxis()
axmy.hist(y_indep, bins=BINS, density=True, orientation="horizontal", color=COLOR_hist)
axmy.set_xlabel("$\\rho_{Y}(y)$ (arb.)")
axmy.set_xticks([])
axmy.set_ylabel("$y$")
axmy.set_yticks([])

# X marginals
axmx1.hist(x_indep, bins=BINS, density=True, color=COLOR_hist)
axmx2.hist(x_dep, bins=BINS, density=True, color=COLOR_hist)

axmx1.set_ylabel("$\\rho_{X}(x)$ (arb.)")
axmx1.set_yticks([])

# compute mutual information
MI_indep = MIEstimate.from_XYKMn_with_RNG(
    X = x_indep.reshape(1, N_samples),
    Y = y_indep.reshape(1, N_samples),
    K = K,
    M = 5,
    n_splits = 5,
    RNG = np.random.default_rng(),
    n_samples = N_samples
)
MI_dep = MIEstimate.from_XYKMn_with_RNG(
    X = x_dep.reshape(1, N_samples),
    Y = y_dep.reshape(1, N_samples),
    K = K,
    M = 5,
    n_splits = 5,
    RNG = np.random.default_rng(),
    n_samples = N_samples
)
print(MI_indep)
print(MI_dep)

for ax, mi in (
    (axh1, MI_indep),
    (axh2, MI_dep),
    ):
    ax.set_xlabel("$x$")
    ax.set_xticks([])
    ax.set_title(
        r"$\hat{\text{I}}_\text{KSG} = " + f"{mi.MI:.3f}" + r" \pm " + f"{mi.std:.3f}" + r"$ nats",
        y = -0.15,
        va="top",
        ha="center",
    )

NO.axis("off")
NOO.axis("off")

plt.savefig("mutual_information_schematic.svg", format="svg", transparent=True, bbox_inches="tight")

    


