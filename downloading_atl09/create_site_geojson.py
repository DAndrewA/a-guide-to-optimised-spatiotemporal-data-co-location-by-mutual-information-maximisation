"""Author: Andrew Martin
Creation date: 31/08/2025

Script to create geojson files containing shapes that are circles surrounding the sites being used in the analysis.
"""


import numpy as np
import cartopy.crs as ccrs
import geojson

SITES = {
    "ny-alesund": dict(lat=78.923, lon=11.922), 
    "hyytiala": dict(lat=61.844, lon=24.287), 
    "juelich": dict(lat=50.908, lon=6.413), 
    "munich": dict(lat=48.148, lon=11.573),
}

def circle_from_location(lat: float, lon:float, N: int = 360, R_m:float = 520*1000):
    """Create an array of latitude-longitude coordinates describing a closed circle around the given latitude and longitude.
    """
    coords_in_tangent_plane = np.array([
        [
            easting:= R_m * np.cos(2 * np.pi * i / N),
            northing:=R_m * np.sin(2 * np.pi * i / N)
        ]
        for i in range(N+1)
    ])
    PROJ = ccrs.Orthographic(
        central_latitude = lat,
        central_longitude = lon
    )
    TRANSFORM = ccrs.PlateCarree()

    x_tangent = coords_in_tangent_plane[:,0]
    y_tangent = coords_in_tangent_plane[:,1]

    transformed = TRANSFORM.transform_points(PROJ, x_tangent, y_tangent)

    coords = [
        [lon, lat]
        for lon, lat in zip(transformed[:,0], transformed[:,1])
    ]
    return coords


def coords_to_geojson(coords: list[list[float]]) -> geojson.Feature:
    # ensure the coordinates are closed
    coords.append(coords[0])
    polygon = geojson.Polygon([coords])
    feature = geojson.Feature(geometry = polygon)
    return feature


def main():
    for site, latlon in SITES.items():
        feature = coords_to_geojson(
            coords = circle_from_location(
                lat = latlon["lat"],
                lon = latlon["lon"],
            )
        )
        fname_out = f"circles/{site}.geojson"
        print(site)
        print(feature)
        print(fname_out)
        with open(fname_out, "w") as f:
            geojson.dump(feature, f)
        print("success")
    print("SUCCESS")


if __name__ == "__main__":
    main()

