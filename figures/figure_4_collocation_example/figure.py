"""Author: Andrew Martin
Creation date: 2/12/25

Script to create Figure 4 in my paper, showing the co-lcoation between ICESat-2 and Cloudnet data.
"""

import sys
sys.path.insert(1,"../..")

from atl09_cloudnet.definitions.collocation import (
    DistanceFromLocation, Duration, RadiusDuration,
    RawATL09, RawCloudnet, 
    CollocationCloudnetATL09, ATL09Event, CloudnetEvent
)
from atl09_cloudnet.definitions import vcf

sys.path.insert(1,"../")
from common.colormaps import COLOR_ATL09, COLOR_Cloudnet
from common.handle_sites import SITES, SITE_locations, SITE_print_names, SITE_argument_names


from sat_val_framework import CollocationEventList

import os
import datetime as dt
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from dataclasses import dataclass, asdict
import pickle

from matplotlib.patches import Circle, Rectangle
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import Divider, Size

from matplotlib.patches import PathPatch
from matplotlib.path import Path
from itertools import product

import cmcrameri.cm as cm

FIG_height = 6
YLIM = [0,10_000]

@dataclass
class PlotArgs:
    c: np.ndarray
    marker: str
    ls: str
    lw: float

pa_atl09 = PlotArgs(
    c = COLOR_ATL09,
    marker = "^",
    ls = "-",
    lw = 2
)
pa_cloudnet = PlotArgs(
    c = COLOR_Cloudnet,
    marker="s",
    ls="-",
    lw=2
)

CMAP_cloudmask = ListedColormap(
    [
        cm.devon(0.9),
        cm.devon(0.5),
        cm.devon(0.1) # attenuation
    ]
)

NORM_cloudmask = BoundaryNorm(
    boundaries = [-0.5, 0.5, 1.5, 2.5],
    ncolors = 3
)

COLOR_spatial_block = cm.grayC(0.6)

# define the inner and outer co-location parametrisations
R_outer_km = 150
R_inner_km = 120

site = "ny-alesund"

tau_outer = dt.timedelta(hours = 18)
tau_inner = dt.timedelta(hours = 6)

outer_parameters = RadiusDuration({
    RawATL09: DistanceFromLocation(
        distance_km = R_outer_km * np.sqrt(2.2), # multiply by sqrt(2) to ensure inscribed square to circle always contains R_outer_km of track
        latitude = SITE_locations[site]["lat"],
        longitude = SITE_locations[site]["lon"],
    ),
    RawCloudnet: Duration(duration = tau_outer)
})

inner_parameters = RadiusDuration({
    RawATL09: DistanceFromLocation(
        distance_km = R_inner_km,
        latitude = SITE_locations[site]["lat"],
        longitude = SITE_locations[site]["lon"],
    ),
    RawCloudnet: Duration(duration = tau_inner)
})

outer_parameters, inner_parameters

#TODO: fix this so it is reproducible away from JASMIN
event = CollocationCloudnetATL09(data={
    RawATL09: ATL09Event(
        fpath = "/gws/ssde/j25b/icecaps/eeasm/paper1/sites/ny-alesund/atl09/ATL09_20210701034659_01151201_006_01_subsetted.h5",
        min_separation_km = 18.4
    ),
    RawCloudnet: CloudnetEvent(
        closest_approach_time = dt.datetime(2021,7,1,4,13,21,513461),
        root_dir = "/gws/ssde/j25b/icecaps/eeasm/paper1/sites/ny-alesund/cloudnet",
        site="ny-alesund",
    ),
})



def make_circle(r):
    t = np.arange(0, np.pi * 2.0, 0.01)
    t = t.reshape((len(t), 1))
    x = r * np.cos(t)
    y = r * np.sin(t)
    return np.hstack((x, y))

def centered_square_with_inscribed_circle(square_halfedge: float, radius: float, alpha=0.5) -> PathPatch:
    assert radius < square_halfedge, "Circle will not be inscribed"

    vertices_square = np.array([
        [-square_halfedge, -square_halfedge],
        [square_halfedge, -square_halfedge],
        [square_halfedge, square_halfedge],
        [-square_halfedge, square_halfedge],
    ])
    codes_square = np.array([Path.MOVETO] + [Path.LINETO]*3, dtype=Path.code_type)

    vertices_circle = make_circle(r=radius)[::-1] # reverse direction of circle path to define handedness of edge
    codes_circle = np.ones(len(vertices_circle), dtype=Path.code_type) * Path.LINETO
    codes_circle[0] = Path.MOVETO

    all_vertices = np.concatenate((vertices_square, vertices_circle))
    all_codes = np.concatenate((codes_square, codes_circle))

    path = Path(
        vertices = all_vertices,
        codes = all_codes,
    )
    print("path created")
    patch = PathPatch(
        path=path,
        facecolor = COLOR_spatial_block,
        alpha = alpha,
        ec=None,
    )
    return patch

def plot_spatial_subset_atl09(collocated_data):
    atl09_subsetter = inner_parameters[RawATL09]

    outer_atl09 = collocated_data[RawATL09]
    inner_atl09 = atl09_subsetter.subset(outer_atl09)

    fig = plt.figure(figsize=(FIG_height,FIG_height), layout="constrained")

    limit_km = R_outer_km
    tick_spacing_km = 100


    projection = ccrs.Orthographic(
        central_latitude = atl09_subsetter.latitude,
        central_longitude = atl09_subsetter.longitude
    )
    ax = plt.axes(
        projection = projection
    )

    # plot coastlines and location of Cloudnet observatory
    ax.coastlines("10m")
    ax.add_feature(cfeature.LAND)
    ax.scatter(
        atl09_subsetter.longitude, atl09_subsetter.latitude, transform=ccrs.PlateCarree(), marker = "*", fc="red", s=500, ec="k", lw=1.5, zorder=10
    )


    # plot the profile ground tracks for the data
    for p, pdisp in zip(outer_atl09.data.profile, (1,2,3,)):
        d = outer_atl09.data.sel(profile=p)
        ax.plot(d["longitude"], d["latitude"], transform=ccrs.PlateCarree(), label = None, lw=1, c="k")

    for p, pdisp in zip(inner_atl09.data.profile, (1,2,3,)):
        d = inner_atl09.data.sel(profile=p)
        ax.plot(d["longitude"], d["latitude"], transform=ccrs.PlateCarree(), label=None, lw = 3, c="green")



    # set the limits according to R_outer_km
    limits = [ limit_km * 1000 * p for p in (-1, 1)]
    ax.set_xlim(limits)
    ax.set_ylim(limits)

    tick_km = [i*tick_spacing_km - tick_spacing_km*(limit_km//tick_spacing_km) for i in range(2*(limit_km//tick_spacing_km) + 1)]
    tick_m = [1000*v for v in tick_km]
    ax.set_xticks(tick_m)
    ax.set_xticklabels(tick_km)
    ax.set_xlabel("Easting (km)")
    ax.set_yticks(tick_m)
    ax.set_yticklabels(tick_km)
    ax.set_ylabel("Northing (km)")


    # circle showing the spatial subsetting of the ATL09 data from the RawATL09 data
    ax.add_patch(
        centered_square_with_inscribed_circle(
            square_halfedge = limit_km*1000,
            radius = atl09_subsetter.distance_km*1000
        )
    )
    ax.add_patch(
        Circle(
            xy = (0,0),
            radius = atl09_subsetter.distance_km*1000,
            lw=1, ls="--", ec="k",
            fc="none", #alpha=0
        )
    )

    return fig, ax


def get_fixed_size_figure_axis() -> (plt.Figure, plt.Axes):
    fig = plt.figure(figsize=(FIG_height, FIG_height/2))

    # size in inches, [padding, axis, padding]
    h = [Size.Fixed(0.55), Size.Fixed(FIG_height - 0.8), Size.Fixed(0.25)]
    v = [Size.Fixed(0.5), Size.Fixed(FIG_height/2 - 0.8), Size.Fixed(0.3)]
    divider = Divider(fig, (0,0,1,1), h, v, aspect=False)

    ax = fig.add_axes(
        divider.get_position(),
        axes_locator=divider.new_locator(nx=1, ny=1)
    )
    return fig, ax


def plot_atl09_feature_mask_data(collocated_data):
    atl09_subsetter = inner_parameters[RawATL09]
    _profile = 1

    outer_atl09 = collocated_data[RawATL09]
    inner_atl09 = atl09_subsetter.subset(outer_atl09)

    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(FIG_height,FIG_height/2), layout="constrained")
    fig, ax = get_fixed_size_figure_axis()

    limit_km = R_outer_km
    limit_subset_km = R_inner_km

    min_time = inner_atl09.data.time.sel(profile=_profile).min().data
    max_time = inner_atl09.data.time.sel(profile=_profile).max().data

    BOUND_lower = outer_atl09.data.time.sel(profile=_profile).min().data
    BOUND_upper = outer_atl09.data.time.sel(profile=_profile).max().data


    cloud_and_attenuation_mask = (outer_atl09.data.sel(profile=_profile).qc_feature_mask % 2 == 1).copy().astype(int)
    cloud_and_attenuation_mask += 2*((outer_atl09.data.sel(profile=_profile).qc_feature_mask // 4) % 2 == 1)
    (
        cloud_and_attenuation_mask
            .plot(
                x="time", y="height",
                ax = ax,
                cmap=CMAP_cloudmask,
                norm=NORM_cloudmask,
                add_colorbar = False,
            )
    )
    # overplot hatched rectangles showing data removed after subsetting
    ax.add_patch(
        Rectangle(
            (BOUND_lower, 0),
            (min_time - BOUND_lower),
            YLIM[1] - YLIM[0],
            fill=False,
            hatch="/",
            ls=""
        )
    )
    ax.add_patch(
        Rectangle(
            (max_time, 0),
            (BOUND_upper - max_time),
            YLIM[1] - YLIM[0],
            fill=False,
            hatch="\\",
            ls=""
        )
    )
    ax.set_ylim(YLIM)
    ax.set_yticks(
        [0, 2000, 4000, 6000, 8000, 10_000],
        ["0", "2", "4", "6", "8", "10"]
    )
    ax.set_ylabel("$z$ (km)")

    ax.axvline(
        min_time,
        ls="--", c="k"
    )
    ax.axvline(
        max_time,
        ls="--", c="k"
    )
    ax.set_title(None)
    ax.set_xlim([BOUND_lower, BOUND_upper])
    ax.set_xlabel("Elapsed GPS seconds, $t$", ha="right")
    ax.set_title("ATL09")
    return fig, ax

def plot_atl09_collocation_criteria_data(collocated_data):
    atl09_subsetter = inner_parameters[RawATL09]
    _profile = 1

    outer_atl09 = collocated_data[RawATL09]
    inner_atl09 = atl09_subsetter.subset(outer_atl09)

    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(FIG_height,FIG_height/2), layout="constrained")
    fig, ax = get_fixed_size_figure_axis()

    limit_km = R_outer_km
    limit_subset_km = R_inner_km

    min_time = inner_atl09.data.time.sel(profile=_profile).min().data
    max_time = inner_atl09.data.time.sel(profile=_profile).max().data

    BOUND_lower = outer_atl09.data.time.sel(profile=_profile).min().data
    BOUND_upper = outer_atl09.data.time.sel(profile=_profile).max().data

    d2s = outer_parameters[RawATL09].get_distance_to_location(outer_atl09).unstack()
    (
        d2s
            .sel(profile=_profile)
            .plot(
                x="time",
                ax = ax,
                **asdict(pa_atl09),
                markevery=(25,75), # a marker every n/25=3 seconds, displaced by 1 second
                markersize=8
            )
    )
    ax.set_ylim([0,None])
    rect_height = (lambda a: a[1] - a[0])( ax.get_ylim() )
    ax.add_patch(
        Rectangle(
            (BOUND_lower, 0),
            (min_time - BOUND_lower),
            rect_height,
            fill=False,
            hatch="/",
            ls=""
        )
    )
    ax.add_patch(
        Rectangle(
            (max_time, 0),
            (BOUND_upper - max_time),
            rect_height,
            fill=False,
            hatch="\\",
            ls=""
        )
    )
    ax.set_ylabel(r"$r\,(t)$ (km)")

    ax.axhline(R_inner_km, ls="--", c="k")
    ax.set_yticks(
        list(np.arange(0,250,50)) + [R_inner_km],
        [str(v) for v in np.arange(0,250,50)] + [r"$R$"]
    )

    # plot vertical lines for the spatial subsetting by co-location criteria
    ax.axvline(
        min_time,
        ls="--", c="k"
    )
    ax.axvline(
        max_time,
        ls="--", c="k"
    )
    ax.set_title(None)
    ax.set_xlim([BOUND_lower, BOUND_upper])
    ax.set_xlabel("Elapsed GPS seconds, $t$", ha="right")
    return fig, ax

def plot_cloudnet_feature_mask_data(collocated_data):
    cloudnet_subsetter = inner_parameters[RawCloudnet]
    
    outer_cloudnet = collocated_data[RawCloudnet]
    inner_cloudnet = cloudnet_subsetter.subset(outer_cloudnet)
    
    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(FIG_height,FIG_height/2), layout="constrained")
    fig, ax = get_fixed_size_figure_axis()
    
    limit_km = tau_outer
    limit_subset_km = tau_inner
    
    min_time = inner_cloudnet.data.time.min().data
    max_time = inner_cloudnet.data.time.max().data
    
    BOUND_lower = outer_cloudnet.data.time.min().data
    BOUND_upper = outer_cloudnet.data.time.max().data
    
    
    (
        outer_cloudnet.data.qc_cloudmask
            .plot(
                x="time", y="height",
                ax = ax,
                cmap=CMAP_cloudmask,
                norm=NORM_cloudmask,
                add_colorbar = False,
            )
    )
    ax.add_patch(
        Rectangle(
            (BOUND_lower, 0),
            (min_time - BOUND_lower),
            YLIM[1] - YLIM[0],
            fill=False,
            hatch="/",
            ls=""
        )
    )
    ax.add_patch(
        Rectangle(
            (max_time, 0),
            (BOUND_upper - max_time),
            YLIM[1] - YLIM[0],
            fill=False,
            hatch="\\",
            ls=""
        )
    )
    ax.set_ylim(YLIM)
    ax.set_yticks(
        [0, 2000, 4000, 6000, 8000, 10_000],
        ["0", "2", "4", "6", "8", "10"]
    )
    ax.set_ylabel("$z$ (km)")
    
    ax.axvline(
        min_time, 
        ls="--", c="k"
    )
    ax.axvline(
        max_time, 
        ls="--", c="k"
    )
    ax.set_xlim([BOUND_lower, BOUND_upper])
    ax.set_xlabel("Time UTC", ha="right")

    ax.set_xticks(
        ax.get_xticks(),
    )
    ax.set_title("Cloudnet")
    
    return fig, ax

def plot_cloudnet_collocation_criteria_data(collocated_data):
    cloudnet_subsetter = inner_parameters[RawCloudnet]
    
    outer_cloudnet = collocated_data[RawCloudnet]
    inner_cloudnet = cloudnet_subsetter.subset(outer_cloudnet)
    
    #fig, ax = plt.subplots(1,1, sharex=True, figsize=(FIG_height,FIG_height/2), layout="constrained")
    fig, ax = get_fixed_size_figure_axis()

    limit_km = tau_outer
    limit_subset_km = tau_inner
    
    min_time = inner_cloudnet.data.time.min().data
    max_time = inner_cloudnet.data.time.max().data
    
    BOUND_lower = outer_cloudnet.data.time.min().data
    BOUND_upper = outer_cloudnet.data.time.max().data
    
    dtau_s = np.abs( outer_cloudnet.data.time - np.datetime64(event[RawCloudnet].closest_approach_time) ) / 1e9
    
    (
        dtau_s
            .plot(
                x="time",
                ax = ax,
                **asdict(pa_cloudnet),
                markevery = (int(dtau_s.count())%180//2,180),
                markersize=8,
            )
    )
    ax.set_ylim([0,None])
    rect_height = (lambda a: a[1] - a[0])( ax.get_ylim() )
    ax.add_patch(
        Rectangle(
            (BOUND_lower, 0),
            (min_time - BOUND_lower),
            rect_height,
            fill=False,
            hatch="/",
            ls=""
        )
    )
    ax.add_patch(
        Rectangle(
            (max_time, 0),
            (BOUND_upper - max_time),
            rect_height,
            fill=False,
            hatch="\\",
            ls=""
        )
    )
    ax.axhline(tau_inner.total_seconds()/2, ls="--", c="k")
    ax.set_ylabel(r"$\left| t - t_{0} \right|$ (hours)")
    ax.set_yticks(
        list(np.arange(0,9,2) * 3600) + [tau_inner.total_seconds()/2],
        [str(v) for v in np.arange(0,9,2)] + [r"$\frac{\tau}{2}$"]
    )
    
    # plot vertical lines for the spatial subsetting by co-location criteria
    ax.axvline(
        min_time, 
        ls="--", c="k"
    )
    ax.axvline(
        max_time, 
        ls="--", c="k"
    )
    ax.set_xlabel("Time UTC", ha="right")

    #ax.set_xticks(
    #    ax.get_xticks(),
    #    ax.get_xticklabels(),
    #    rotation=20
    #)
    ax.set_xlim([BOUND_lower, BOUND_upper])
    ax.set_title(None)
    
    return fig, ax

def plot_homogenised_data(collocated_data):
    homgenised_data = collocated_data.homogenise_to(vcf.VCF_240m)

    fig, ax = plt.subplots(1,1, figsize=(FIG_height/2, FIG_height), layout="constrained")

    homgenised_data[RawATL09].data.plot(
        y="height", ax=ax,
        **asdict(pa_atl09),
        markevery=2,
        markersize=8,
        label="ATL09",
    )
    homgenised_data[RawCloudnet].data.plot(
        y="height", ax=ax,
        **asdict(pa_cloudnet),
        label="Cloudnet",
        markevery=(1,2),
        markersize=8,
    )

    ax.set_ylabel("$z$ (km)")
    ax.set_xlabel(r"$\text{VCF}$")
    ax.set_xlim([-0.01,1])
    ax.set_ylim(YLIM)
    ax.set_yticks(
        [0, 2000, 4000, 6000, 8000, 10_000],
        ["0", "2", "4", "6", "8", "10"]
    )

    ax.set_xticks([0,0.5,1])

    ax.legend()
    return fig, ax


print(f"loading data")
collocated_data = event.load_with_joint_parameters(outer_parameters)
print("data loaded succesfully")
# svg-able
for plot_func, savename in (
#    (plot_spatial_subset_atl09,"spatial_subset"), 
    (plot_atl09_collocation_criteria_data, "atl09_criteria"), 
    (plot_cloudnet_collocation_criteria_data, "cloudnet_criteria"), 
#    (plot_homogenised_data, "VCF"),
):
    f,a = plot_func(collocated_data)
    f.savefig(f"{savename}.svg", format="svg", transparent=True)

# png
for plot_func, savename in (
    (plot_atl09_feature_mask_data, "atl09_feature"),
    (plot_cloudnet_feature_mask_data, "cloudnet_feature"),
):
    f,a = plot_func(collocated_data)
    f.savefig(f"{savename}.png", dpi=600, transparent=True)


def plot_colorbar():
    width_reduction_factor = 0.68
    aspect = 20

    fig, cax = plt.subplots(
        1,
        1,
        figsize=(FIG_height/2, FIG_height/2),
        layout="constrained",
        #gridspec_kw = dict(right=0, left=-50)
    )

    mappable = ScalarMappable(cmap=CMAP_cloudmask, norm=NORM_cloudmask)
    mappable.set_array([])
    plt.colorbar(mappable, cax=cax)
    cax.set_yticks(
        [0, 1, 2],
        ["no cloud", "cloud", "attenuation"],
        rotation=-90,
        va="center",
        #va = ["center", "top", "top"]
        ha="left"
    )
    cax.yaxis.set_ticks_position('right')
    cax.set_box_aspect(aspect)
    return fig, cax

f,a = plot_colorbar()
f.savefig(
    "colorbar.svg",
    format="svg",
    transparent=True,
    bbox_inches="tight"
)


