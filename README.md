# Example: ATL09 and Cloudnet data

The code in this repository will demonstrate how the framework can be used by applying it to a validation of ICESat-2 ATL09 data against Cloudnet data.
The ATL09 data will be collocated with the Cloudnet data by taking vertical profiles that fall within a radius R of the Cloudnet site, and Cloudnet data within a window of duration tau, centered on the time of closest approach.

In order to run the analysis, we implement:
+ `RawATL09 (RawData)` to handle loading ATL09 data from (possibly multiple) `.h5` files.
+ `RawCloudnet (RawData)` to handle loading Cloudnet data from (possibly multiple) `.nc` files.
+ `RadiusDuration (CollocationParametrisation)` to handle passing collocation parameters to subset the raw data.
+ `VCFProfile (HomogenisedData)` to create homogenised Vertical Cloud Fraction profiles from raw data.
+ `CollocationCloudnetATL09 (CollocationEvent)` to handle all of the information required to load raw data for a collocation event.
+ `SchemeCloudnetATL09RadiusDuration (CollocationScheme)` to handle generating `CollocationEvent`s using definitions for `RawData`.
