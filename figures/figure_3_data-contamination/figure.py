"""Author: Andrew Martin
Creation date: 2/12/25

Script to handle creating panels for figure showing the change from being data limited to having samples contaminated with independent data.
"""

import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(1, "../")
from common.colormaps import CMAP_probability, MAPPABLE_probability
from common.synthetic_data import (
    gen_xy_samples_mixture,
    x_of_u_given_F_of_x,
    F_a_of_x,
    F_b_of_x,
)
from common.handle_MI_datasets import K

sys.path.insert(1, "../../")
from atl09_cloudnet.holmes_et_al_2019.python_interface import MIEstimate


base = np.exp(1)
RHOSTAR = np.sqrt(1 - np.power(base, -2)) # almost guarantees MI = 1 nat
RHOSTAR = 0.99

BINS = np.linspace(0,1,101)

N_star = 3000
kappa_star = 1
kappa_of_r = lambda r: kappa_star if r < 1 else kappa_star/r**2
n_of_r = lambda r: int(N_star*r**2)
R = (0.2,1,2,4)

# generate one-to-one mapping of X to Y
one_to_one_u = np.linspace(0.01, 0.99, 500)
(_, one_to_one_x) = x_of_u_given_F_of_x(u=one_to_one_u, F_of_x=F_a_of_x)
(_, one_to_one_y) = x_of_u_given_F_of_x(u=one_to_one_u, F_of_x=F_b_of_x)
one_to_one_y = 1 - one_to_one_y



fig, axs = plt.subplots(1,5, figsize=(13,3), width_ratios=[3,3,3,3,1])

for ax, r in zip(axs, R):
    mixture = kappa_of_r(r)
    n_samples = n_of_r(r)

    (u,v), (x,y) = gen_xy_samples_mixture(
        n_samples = n_samples,
        rhostar = RHOSTAR,
        mixture = mixture,
    )
    y = 1 - y

    ax.hist2d(x,y, density=True, cmap=CMAP_probability, bins=BINS)

    ax.plot(one_to_one_x, one_to_one_y, ls="--", c="w", lw=1, alpha=0.5)

    mi = MIEstimate.from_XYKMn_with_RNG(
        X = x.reshape(1, n_samples),
        Y = y.reshape(1, n_samples),
        K = K,
        M = 5,
        n_splits = 5,
        RNG = np.random.default_rng(),
        n_samples = n_samples,
    )

    ax.set_xlabel("$x$")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_title(
        r"$\hat{\text{I}}_\text{KSG}=" + f"{mi.MI:.3f}" + r"\pm" + f"{mi.std:.3f}" + r"$ nats",
        y=-0.15,
        va="top",
        ha="center"
    )
    ax.set_box_aspect(1)
    ax.set_xlim([0,1])
    ax.set_xlim([0,1])
    print(f"{mi.MI:.3f} pm {mi.std:.3f}")
axs[0].set_ylabel("$y$")

cax = axs[-1]
plt.colorbar(MAPPABLE_probability, cax=cax, fraction=0.7)
cax.set_ylabel("Normalised probability density")
cax.set_yticks([])
cax.yaxis.set_ticks_position("left")
cax.set_box_aspect(12)

plt.savefig("mutual_information_contamination.svg", format="svg", transparent=True)


