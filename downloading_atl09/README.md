# Downloading ATL09 data

The scripts in this folder describe how the ATL09 data for the analysis was obtained.

Firstly, `geojson` files are created as circles centered on the Cloudnet sites, as defined by `create_site_geojson.py`.

Then, the ![NASA Harmony API](https://harmony.earthdata.nasa.gov/docs) is used to construct HTTPS POST requests to subset ATL09 data based on the relevant shapefiles.

# obtaining data
1. Run `python create_site_geojson.py` to create geojson files describing circles around the Cloudnet sites.
2. Run `python generate_harmony_requests.py` and select option `Y` for all sites to generate requests for.
3. Output the job ids displayed from `generate_harmony_requests.py` to the file `harmony-job-ids`
4. Run `. urls_by_id/get_job_urls` to obtain urls for the files to download.
5. On the xfer servers, and in the target directory `$MI_MAXIMISATION_RAW_DATA_DIRECTORY/sites/<site>/atl09`, run the command `nargs -n 1 curl -O -b ~/.urs_cookies -c ~/.urs_cookies -L -n < <urls-file>` to download the ATL09 data (recommended to run in tmux as this process will be long-running)
