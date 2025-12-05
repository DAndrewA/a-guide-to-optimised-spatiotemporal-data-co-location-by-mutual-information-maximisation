"""Author: Andrew Martin
Creation date: 5/12/25

Script to generate the Figure 12 demonstration of normalised across-track orbital densities
"""

import sys
sys.path.insert(1,"../")
from common.rho_orbits import (
    normalised_orbital_density_degrees
)

import cmcrameri.cm as cm
import numpy as np

import matplotlib.pyplot as plt

INCLINATIONS = {
    "ISS": 51.64,
    "EarthCARE": 97.05,
    "A-TRAIN": 98.14,
    "ICESat-2": 92.00,
}

COLORS = {
    k: cm.batlowS(2*i)
    for i, k in enumerate(INCLINATIONS.keys())
}

acute_angle = lambda theta: abs(theta) if abs(theta) <= np.pi/2 else np.pi - np.abs(theta)
acute_angle_degrees = lambda theta: np.rad2deg( acute_angle(np.deg2rad(theta)) )

f, a = plt.subplots(1,1, figsize=(8.3,4))

for satellite, inclination in INCLINATIONS.items():
    color = COLORS[satellite]
    max_lat = acute_angle_degrees(inclination)
    print(satellite, max_lat)

    latitudes = np.linspace(0, max_lat, 5000)

    slice_all_but_final = slice(None, -1)
    #print(np.min(latitudes), np.max(latitudes), latitudes.shape)
    a.plot(
        latitudes[slice_all_but_final],
        normalised_orbital_density_degrees(
            latitude_1 = latitudes,
            latitude_2 = max_lat,
        )[slice_all_but_final],
        c=color,
        label=satellite,
        marker = "o" if satellite=="ICESat-2" else None,
        markevery=(50,250)
    )
    a.axvline(max_lat, ls="--", c=color)
a.legend()
a.set_ylim([0,15])
a.set_xlim([0,90])
for rho in (1,2,5,10): a.axhline(rho, ls=(-3*rho,(15,20)), c="k", lw=0.5)
a.set_yticks([0,2,5,10,15])
a.set_xlabel(r"latitude ($\degree$N)")
a.set_ylabel(r"$\rho_{\text{orbits}}$")

plt.savefig("rho_orbits.svg", format="svg", transparent=True, bbox_inches="tight")

