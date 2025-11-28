"""Author: Andrew Martin
Creation date: 28/11/25

Script to produce 3 figures for the mutual information bounds appendix
"""

import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, "../../")
from atl09_cloudnet.holmes_et_al_2019.python_interface import MIEstimate

import scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects

K=10


# rhostar chosen to maintain mutual information = 1 in units of given base
base = np.exp(1)
RHOSTAR = np.sqrt(1 - np.power(base,-2))#np.sqrt(3)/2

FIGSIZE = (5,3)

def gen_samples_mixture(n_samples: int, rhostar: float, mixture: float, rng: np.random.Generator | None = None):
    assert rhostar >= -1 and rhostar <= 1
    assert 0 <= mixture and mixture <= 1
    assert isinstance(n_samples, int)
    if rng is None:
        rng = np.random.default_rng()

    if mixture == 0: return rng.uniform(low=0, high=1, size=(n_samples, 2))

    n_gaussian = int(np.ceil(n_samples * mixture))

    cov_gaussian = np.array([
        [1, rhostar],
        [rhostar, 1]
    ])
    xy_gaussian = scipy.stats.multivariate_normal.rvs(cov=cov_gaussian, size=n_gaussian)
    uv_gaussian = 0.5*( 1 + scipy.special.erf(xy_gaussian / np.sqrt(2)) )

    if mixture == 1: return uv_gaussian

    n_independent = n_samples - n_gaussian
    uv_independent = rng.uniform(low=0, high=1, size=(n_independent, 2))

    uv = np.concat([uv_gaussian, uv_independent])
    return uv


def figure_1():
    N_values = np.logspace(6,2,20).astype(int)
    print(N_values)
    MI_estimates_dependent = [
        (print(N),
        MIEstimate.from_XYKMn_with_RNG(
            X = (
                uv:=gen_samples_mixture(
                    n_samples = int(N),
                    rhostar=RHOSTAR,
                    mixture=1,
                )
            )[:,0].reshape(1,N),
            Y = uv[:,1].reshape(1,N),
            K = K,
            M = 5,
            n_splits=5,
            RNG = np.random.default_rng(),
            n_samples=int(N)
        ))[-1]
        for N in N_values
    ]
    
    MI_estimates_independent = [
        (print(N),
        MIEstimate.from_XYKMn_with_RNG(
            X = (
                uv:=gen_samples_mixture(
                    n_samples = int(N),
                    rhostar=0,
                    mixture=0,
                )
            )[:,0].reshape(1,N),
            Y = uv[:,1].reshape(1,N),
            K = K,
            M = 5,
            n_splits=5,
            RNG = np.random.default_rng(),
            n_samples=int(N)
        ))[-1]
        for N in N_values
    ]

    fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")

    plotN = 1/N_values#-np.log10(N_values)

    ax.axhline(0, ls="--", c="k", lw=1, alpha=0.5)
    ax.axhline(MI_estimates_dependent[0].MI, ls="--", c="k", lw=1, alpha=0.5)

    ax.errorbar(
        x=plotN,
        y=[
            MIE.MI
            for MIE in MI_estimates_dependent
        ],
        yerr=[
            MIE.std
            for MIE in MI_estimates_dependent
        ],
        ls="none",
        label="dependent",
        marker="x",
        capsize=2,
    )

    ax.errorbar(
        x=plotN,
        y=[
            MIE.MI
            for MIE in MI_estimates_independent
        ],
        yerr=[
            MIE.std
            for MIE in MI_estimates_independent
        ],
        ls="none",
        label="independent",
        marker="x",
        capsize=2,
    )

    ax.set_xscale("log")
    ax.set_xlabel(r"$1 \slash N$")
    ax.legend()
    ax.set_ylabel(r"$\hat{\text{I}}_\text{KSG}$ (nats)")
    
    plt.savefig("mi_as_function_of_n.svg", transparent=True, format="svg", bbox_inches="tight")

def figure_2():
    N = 100_000

    mixing_ratios = np.linspace(0,1,21)

    MI_estimates_mixtures = [
        (print(mixing_ratio),
        MIEstimate.from_XYKMn_with_RNG(
            X = (
                uv:=gen_samples_mixture(
                    n_samples = int(N),
                    rhostar=RHOSTAR,
                    mixture=mixing_ratio,
                )
            )[:,0].reshape(1,N),
            Y = uv[:,1].reshape(1,N),
            K = K,
            M = 5,
            n_splits=5,
            RNG = np.random.default_rng(),
            n_samples=int(N)
        ))[-1]
        for mixing_ratio in mixing_ratios
    ]

    fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")

    ax.plot(
        [0,1],[1,0], 
        ls=(0,(10,2)), c="k", lw=1, 
        label="theoretical bound",
        path_effects = [patheffects.withTickedStroke(spacing=12, angle=-90, length=0.3)],
    )
    ax.axhline(0, ls="--", c="k", lw=1, alpha=0.5)

    ax.errorbar(
        x=1 - mixing_ratios,
        y=[
            MIE.MI
            for MIE in MI_estimates_mixtures
        ],
        yerr=[
            MIE.std
            for MIE in MI_estimates_mixtures
        ],
        label=r"$\hat{\text{I}}_\text{KSG} ( X;Y \, | \, \kappa )$",
        marker="x",
        capsize=2,
    )

    ax.set_xlim([0,1])
    ax.set_ylim([None, 1])

    ax.set_xticks([0,1], [1,0])
    ax.set_yticks([0,1])

    ax.set_xlabel(r"$\kappa$")
    ax.set_ylabel(r"$\hat{\text{I}}_\text{KSG}$ (nats)")

    ax.legend()
    ax.set_box_aspect(1)

    plt.savefig("mi_as_function_of_mixing_ratio.svg", format="svg", transparent=True, bbox_inches="tight")

def figure_3():
    R_over_Rstar_values = np.concat([np.arange(0.2,1,0.2), np.linspace(1,4,11)])

    kappa_of_R = lambda R: np.ones_like(R)*(R <= 1) + np.ones_like(R)/np.power(R,2)*(R>1)
    N_star = [
        100,
        1_000,
        10_000,
        100_000,
    ]
    sampling_densities = [
        Nstar / np.pi
        for Nstar in N_star
    ]
    N_of_R = lambda R, sampling_density: np.pi * sampling_density * np.power(R,2)

    def safe_MI(**kwargs):
        try:
            return MIEstimate.from_XYKMn_with_RNG(**kwargs)
        except:
            return MIEstimate(
                MI=0,
                std=0,
                N=N,
                sigma_i=list(),
                M=0,
                n_splits=0,
                K=K
            )

    MI_by_R_by_sampling_density = [
        [
            (
                N:=int(N_of_R(R,sampling_density)),
                mixture:=kappa_of_R(R),
                uv:=np.atleast_2d(gen_samples_mixture(
                    n_samples=N,
                    rhostar=RHOSTAR,
                    mixture=mixture,
                )),
                print(sampling_density, R, N, mixture, uv.shape),
                safe_MI(
                    X = uv[:,0].reshape(1,N),
                    Y = uv[:,1].reshape(1,N),
                    K = K,
                    M = 5,
                    n_splits=5,
                    RNG = np.random.default_rng(),
                    n_samples=N
                )
            )[-1]
            for R in R_over_Rstar_values
        ]
        for sampling_density in sampling_densities
    ]

    fig, ax = plt.subplots(1,1, figsize=FIGSIZE, layout="constrained")

    theoretical_bound = kappa_of_R(R_over_Rstar_values)
    markers = ["^","s","o","*"]

    theoretical_line = ax.plot(
        [0,*R_over_Rstar_values],
        [1,*theoretical_bound],
        ls=(0,(10,2)), c="k", lw=1,
        path_effects = [patheffects.withTickedStroke(spacing=12, angle=-90, length=0.3)],
        label="theoretical bound"
    )

    ax.axhline(0, ls="--", c="k", lw=1, alpha=0.5)
    ax.axvline(1, ls="--", c="k", lw=1, alpha=0.5)
    #ax.axhline(1, ls="--", c="k", lw=1, alpha=0.5)

    for i, Nstar in enumerate(N_star):
        MIE_list = [ MIE for MIE in MI_by_R_by_sampling_density[i] ]

        ax.errorbar(
            x=R_over_Rstar_values,
            y=[
                MIE.MI
                for MIE in MIE_list
            ],
            yerr=[
                MIE.std
                for MIE in MIE_list
            ],
            label=f"$N^* = 10^{np.log10(Nstar):.0f}$",
            marker = markers[i],
            markevery = (i,4),
            markersize=8,
            capsize=2,
        )

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(r"$\hat{\text{I}}_\text{KSG} (X;Y \, | \, r<R)$ (nats)")

    ax.set_xticks(
        [0,1,2,3,4],
        [0, r"$R^*$",r"$2R^*$",r"$3R^*$",r"$4R^*$"]
    )
    ax.set_yticks([0,1])
    ax.set_xlim([0,4])
    ax.set_ylim([None,None])

    ax.legend()
    plt.savefig("mi_as_function_of_R.svg", format="svg", transparent=True, bbox_inches="tight")


figure_1()
figure_2()
figure_3()
