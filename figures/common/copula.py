"""Author: Andrew Martin
Creation date: 27/11/25

Functions pertaining to the generation and use of bivariate copulas
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects

from .common import get_figure_panel
from .colormaps import CMAP_copula, NORM_copula

from typing import Self
from dataclasses import dataclass


COPULA_bins = 25


def fill_non_finite(da: xr.DataArray, v=np.nan) -> xr.DataArray:
    return da.copy().where(np.isfinite(da), v)


@dataclass
class BivariateCopula:
    data: xr.DataArray

    @classmethod
    def generate(cls, data_X: np.ndarray, data_Y: np.ndarray, nbins_nominal: int=COPULA_bins) -> Self:
        """Generates a bi-variate copula for data data_X and data_Y.
        The number of bins supplied is nominal, as degenerate values in the data lead to discontinuities in the cumulative distribution, inducing delta distributions in the copula density.
        Thus, degeneracies in the inverse-mapping x = F^{-1}(u) are lifted, by generating copula with u and v coordinates only at the maximum value of the degenerate data.
        """
        desired_u = np.linspace(0,1, nbins_nominal + 1)
        desired_v = np.linspace(0,1, nbins_nominal + 1)

        quantiles_X = np.quantile(data_X, desired_u)
        quantiles_Y = np.quantile(data_Y, desired_v)

        # determine the x and y quantiles which can be used to lift the degeneracy of the distribution at the lower end
        uqX, iuqX = np.unique(quantiles_X, return_index=True)
        # the unique indices, iuqX, are the lowest index of a given unique value. We want the highest. Hence, iuqX[1:]-1
        u_non_degen = np.concat((
            [0],
            desired_u[iuqX[1:]-1],
            desired_u[-1:],
        ))

        uqY, iuqY = np.unique(quantiles_Y, return_index=True)
        v_non_degen = np.concat((
            [0],
            desired_v[iuqY[1:]-1],
            desired_v[-1:],
        ))

        # generate the data that defines the copula
        V,U = np.meshgrid(v_non_degen, u_non_degen)
        Y,X = np.meshgrid(
            np.insert(uqY,[0], [uqY.min()-1]),
            np.insert(uqX,[0], [uqX.min()-1]),
        )

        d = np.logical_and(
            (np.expand_dims(data_X,(-1,-2)) <= X),
            (np.expand_dims(data_Y,(-1,-2)) <= Y)
        ).mean(axis=0)

        copula = xr.DataArray(
            data = d,
            #coords = {
            #    "u": bins_copula,
            #    "v": bins_copula
            #},
            coords = {
                "u": u_non_degen,
                "v": v_non_degen,
            },
            dims = ["u", "v"]
        ).rename("copula")
        return cls(data=copula)


    @property
    def density(self) -> xr.DataArray:
        return (
            fill_non_finite(
                self.data
                    .differentiate(coord="u")
                    .differentiate(coord="v")
            )
            .clip(min=0)
        )


    @property
    def RMSD(self) -> float:
        return float(
            np.sqrt(
                np.power(self.density-1, 2).mean(skipna=True)
            )
        )


    @property
    def cmin(self) -> float:
        return float(self.density.min(skipna=True))


    @property
    def cmax(self) -> float:
        return float(self.density.max(skipna=True))


    @property
    def c11(self) -> float:
        return float(self.density.sel(dict(u=1, v=1)))


    def plot_density(self, cbar: bool=False) -> (plt.Figure, plt.Axes):
        fig, ax = get_figure_panel()

        density = self.density
        density.plot(
            x = "u",
            y = "v",
            cmap = CMAP_copula,
            norm = NORM_copula,
            ax = ax,
            add_colorbar = cbar,
        )
        contour_copula = density.plot.contour(
            x = "u",
            y = "v",
            levels = [1],
            colors = "k",
            ax = ax,
        )
        contour_copula.set(path_effects = [
            patheffects.withTickedStroke(length=0.5, spacing=5, angle=90)
        ])

        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.set_box_aspect(1)

        return fig, ax
