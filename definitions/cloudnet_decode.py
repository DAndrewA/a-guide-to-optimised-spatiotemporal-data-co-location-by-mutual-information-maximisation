"""
Created on Monday 25 Nov 2024

@author: heather guy

Functions for parsing and plotting cloudnet
categorization and classification data.
"""

import xarray as xr
import numpy as np
from matplotlib.colors import ListedColormap

# Plotting colormaps
_COLORS = {
    "green": "#3cb371",
    "darkgreen": "#253A24",
    "lightgreen": "#70EB5D",
    "yellowgreen": "#C7FA3A",
    "yellow": "#FFE744",
    "orange": "#ffa500",
    "pink": "#B43757",
    "red": "#F57150",
    "shockred": "#E64A23",
    "seaweed": "#646F5E",
    "seaweed_roll": "#748269",
    "white": "#ffffff",
    "lightblue": "#6CFFEC",
    "blue": "#209FF3",
    "skyblue": "#CDF5F6",
    "darksky": "#76A9AB",
    "darkpurple": "#464AB9",
    "lightpurple": "#6A5ACD",
    "purple": "#BF9AFF",
    "darkgray": "#2f4f4f",
    "lightgray": "#ECECEC",
    "gray": "#d3d3d3",
    "lightbrown": "#CEBC89",
    "lightsteel": "#a0b0bb",
    "steelblue": "#4682b4",
    "mask": "#C8C8C8",
}
class_colors=(
        ("_Clear sky", _COLORS["white"]),
        ("Droplets", _COLORS["lightblue"]),
        ("Drizzle or rain", _COLORS["blue"]),
        ("Drizzle & droplets", _COLORS["purple"]),
        ("Ice", _COLORS["lightsteel"]),
        ("Ice & droplets", _COLORS["darkpurple"]),
        ("Melting ice", _COLORS["orange"]),
        ("Melting & droplets", _COLORS["yellowgreen"]),
        ("Aerosols", _COLORS["lightbrown"]),
        ("Insects", _COLORS["shockred"]),
        ("Aerosols & insects", _COLORS["pink"]),
        ("No data", _COLORS["mask"])
    )

detection_colors=(
        ("_Clear sky", _COLORS["white"]),
        ("Lidar only", _COLORS["yellow"]),
        ("Uncorrected atten.", _COLORS["seaweed_roll"]),
        ("Radar & lidar", _COLORS["green"]),
        ("_No radar but unknown atten.", _COLORS["purple"]),
        ("Radar only", _COLORS["lightgreen"]),
        ("_No radar but known atten.", _COLORS["orange"]),
        ("Corrected atten.", _COLORS["skyblue"]),
        ("Clutter", _COLORS["shockred"]),
        ("_Lidar molecular scattering", _COLORS["pink"]),
        ("No data", _COLORS["mask"]))

# Define categories for classification
class_names= [cc[0] for cc in class_colors]
class_cols=[cc[1] for cc in class_colors]
# Using a predefined colormap (e.g., Set3)
class_cmap = ListedColormap(class_cols[:-1])
class_cmap.set_bad(color=class_cols[-1])

stat_names= [cc[0] for cc in detection_colors]
stat_cols=[cc[1] for cc in detection_colors]
stat_cmap = ListedColormap(stat_cols[:-1])
stat_cmap.set_bad(color=stat_cols[-1])

def decode_cat(arr):
    # Bit 0: Small liquid droplets are present.
    # Bit 1: Falling hydrometeors are present; if Bit 2 is set then these are most
    # likely ice particles, otherwise they are drizzle or rain drops.
    # Bit 2: Wet-bulb temperature is less than 0 degrees C, implying the phase of
    # Bit-1 particles.
    # Bit 3: Melting ice particles are present.
    # Bit 4: Aerosol particles are present and visible to the lidar.
    # Bit 5: Insects are present and visible to the radar.
    arr_flat = arr.flatten()
    nan_mask = np.isnan(arr_flat)
    arr_nonan = arr_flat[~nan_mask]
    arr_decode=[bin(int(c))[2:].zfill(6) for c in arr_nonan]
    liq = np.asarray([int(decode[-1]) for decode in arr_decode ])
    falling = np.asarray([int(decode[-2]) for decode in arr_decode ])
    cold = np.asarray([int(decode[-3]) for decode in arr_decode ])
    melting = np.asarray([int(decode[-4]) for decode in arr_decode ])
    aerosol = np.asarray([int(decode[-5]) for decode in arr_decode ])
    insects = np.asarray([int(decode[-6]) for decode in arr_decode ])

    # reinsert nans and reshape
    liq = replace_nan(liq,nan_mask).reshape(np.shape(arr))
    falling = replace_nan(falling,nan_mask).reshape(np.shape(arr))
    cold = replace_nan(cold,nan_mask).reshape(np.shape(arr))
    melting = replace_nan(melting,nan_mask).reshape(np.shape(arr))
    aerosol = replace_nan(aerosol,nan_mask).reshape(np.shape(arr))
    insects = replace_nan(insects,nan_mask).reshape(np.shape(arr))

    return liq,falling,cold,melting,aerosol,insects


def recode_cat(liq,falling,cold,melting,aerosol,insects):
    # flatten the input arrays
    l=liq.flatten()
    f=falling.flatten()
    c=cold.flatten()
    m=melting.flatten()
    a=aerosol.flatten()
    ins=insects.flatten()

    # Make the binary strings
    bi_str=['%i%i%i%i%i%i'%(ins[i],a[i],m[i],c[i],f[i],l[i]) for i in range(0,len(l)) ]
    # cal the integers
    ints = [int(s,2) for s in bi_str]
    # reshape
    cat_bits=np.asarray(ints).reshape(np.shape(liq))
    return cat_bits

def decode_quality(arr):
    # Bit 0: An echo is detected by the radar.
    # Bit 1: An echo is detected by the lidar.
    # Bit 2: The apparent echo detected by the radar is ground clutter or some
    #        other non-atmospheric artifact.
    # Bit 3: The lidar echo is due to clear-air molecular scattering.
    # Bit 4: Liquid water cloud, rainfall or melting ice below this pixel
    #        will have caused radar and lidar attenuation; if bit 5 is set then
    #        a correction for the radar attenuation has been performed;
    #        otherwise do not trust the absolute values of reflectivity factor.
    #        No correction is performed for lidar attenuation.
    # Bit 5: Radar reflectivity has been corrected for liquid-water attenuation
    #        using the microwave radiometer measurements of liquid water path
    #        and the lidar estimation of the location of liquid water cloud;
    #        be aware that errors in reflectivity may result.
    arr_flat = arr.flatten()
    nan_mask = np.isnan(arr_flat)
    arr_nonan = arr_flat[~nan_mask]
    arr_decode=[bin(int(c))[2:].zfill(6) for c in arr_nonan]
    radar = np.asarray([int(decode[-1]) for decode in arr_decode ])
    lidar = np.asarray([int(decode[-2]) for decode in arr_decode ])
    radar_clutter = np.asarray([int(decode[-3]) for decode in arr_decode ])
    lidar_clearair = np.asarray([int(decode[-4]) for decode in arr_decode ])
    attenuated = np.asarray([int(decode[-5]) for decode in arr_decode ])
    atten_corrected = np.asarray([int(decode[-6]) for decode in arr_decode ])

    # reinsert nans and reshape
    radar = replace_nan(radar,nan_mask).reshape(np.shape(arr))
    lidar = replace_nan(lidar,nan_mask).reshape(np.shape(arr))
    radar_clutter = replace_nan(radar_clutter,nan_mask).reshape(np.shape(arr))
    lidar_clearair = replace_nan(lidar_clearair,nan_mask).reshape(np.shape(arr))
    attenuated = replace_nan(attenuated,nan_mask).reshape(np.shape(arr))
    atten_corrected = replace_nan(atten_corrected,nan_mask).reshape(np.shape(arr))

    return radar,lidar,radar_clutter,lidar_clearair,attenuated,atten_corrected


def replace_nan(nonan,mask):
    #if non nan is a flat array with nans removed
    # and mask is the orginal masked array where nans are
    # add the nans back into nonan and return.
    c =  np.empty_like(mask)
    c.fill(np.nan)
    c[~mask] = nonan
    return c


"""NOTE: all code written above was written by Heather Guy. All proceeding code is written by Andrew Martin, for the purposes of interfacing the above functions with xarray.
Date: 22/1/25
"""

def decode_cat_ds(ds_cat: xr.Dataset):
    """Function to decode the category_bits field of a cloudnet categorize data product"""
    ds_decoded = ds_cat.copy()
    liq,falling,cold,melting,aerosol,insects = decode_cat(ds_cat.category_bits.values)
    radar,lidar,radar_clutter,lidar_clearair,attenuated,atten_corrected = decode_quality(ds_cat.quality_bits.values)
    ds_decoded = ds_decoded.assign({
        name: (("time", "height"), data)
        for name, data in zip(
            ("liq","falling","cold","melting","aerosol","insects","radar","lidar","radar_clutter","lidar_clearair","attenuated","atten_corrected"),
            (liq,falling,cold,melting,aerosol,insects,radar,lidar,radar_clutter,lidar_clearair,attenuated,atten_corrected),
        )
    })
    return ds_decoded
