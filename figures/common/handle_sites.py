"""Dictionaries containing the sites, and translation between code-use site labels and latex names/function argument names"""


SITES = ("ny-alesund", "hyytiala", "juelich", "munich",)

SITE_locations = {
    "ny-alesund": dict(lat=78.923, lon=11.922), 
    "hyytiala": dict(lat=61.844, lon=24.287), 
    "juelich": dict(lat=50.908, lon=6.413), 
    "munich": dict(lat=48.148, lon=11.573),
}

SITE_print_names = {
    "ny-alesund": r"Ny-\r{A}lensund",
    "hyytiala": "Hyytiala",
    "juelich": "J\\\"ulich",
    "munich": "Munich",
}

SITE_argument_names = {
    site: site.replace("-","_")
    for site in SITES
}