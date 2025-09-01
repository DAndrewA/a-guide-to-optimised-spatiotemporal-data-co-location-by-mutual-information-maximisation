"""Author: Andrew Martin
Creation date: 31/08/2025

Script to generate requests that can be ran to submit jobs to NASA Harmony.
"""

import urllib.parse
import datetime as dt
from harmony import Request, Collection, WKT, Client
import geojson



TEMPORAL = {
    "start": dt.datetime(2018,10,1),
    "stop": dt.datetime(2025,1,1)
}
# FOR TESTING PURPOSES
#TEMPORAL = {
#    "start": dt.datetime(2020,1,1),
#    "stop": dt.datetime(2020,4,1),
#}
collection_id = "C2649212495-NSIDC_CPRD"



def polygon_coords_to_wkt(coords: list[list[float]]) -> str:
    return "".join([
        "POLYGON((",
        ", ".join([
            f"{lon:.6f} {lat:.6f}"
            for [lon, lat] in coords
        ]),
        "))"
    ])



def geojson_from_site(site: str) -> geojson.Feature:
    with open(f"circles/{site}.geojson", "r") as f:
        feature = geojson.load(f)
    
    return feature



def request_from_wkt(wkt: str):
    request = Request(
        collection = Collection(id=collection_id),
        spatial = WKT(wkt=wkt),
        temporal = TEMPORAL
    )
    return request



SITES = ("ny-alesund", "hyytiala", "juelich", "munich",)


def main():
    requests_by_site = dict()
    for site in SITES:
        feature = geojson_from_site(site)

        wkt = polygon_coords_to_wkt(
                coords = feature.geometry.coordinates[0][:-1] # final closing vertex removed
        )

        request = request_from_wkt(wkt)

        requests_by_site[site] = request

    print(requests_by_site)

    client = Client()
    for site, req in requests_by_site.items():
        if (inp:=input(f"Type Y to submit for {site} | Q to quit | anything else to skip: ")) == "Y":
            job_id = client.submit(req)
            print(f"Job id for {site} is {job_id}")
        elif inp == "Q": break


if __name__ == "__main__":
    main()
