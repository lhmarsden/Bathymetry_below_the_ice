#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import sys

# Traverse THREDDS server, isolate 'N' for Greenland.
url = 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/15s/15s_bed_elev_netcdf/ETOPO_2022_v1_15s_N75W030_bed.nc'

ds = xr.open_dataset(url)

bathymetry = ds['z']

# Subset the data to the area of interest
lat_min = 72
lat_max = 74
lon_min = -17
lon_max = -16
# Elevation range for colour scale
vmax = 3000
vmin = vmax * -1
# Contour interval
contour_interval = 500
contour_levels = np.arange(vmin, vmax + contour_interval, contour_interval)
# Plot only every nth sample in both lat and lon to speed up processing
sampling_factor = 10

bathymetry_zoomed = bathymetry.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))

# Check if there is data in the subset
if bathymetry_zoomed.size == 0:
    print("No data in the specified range for file")
    sys.exit()
else:
    bathymetry = bathymetry_zoomed

# Select every nth sample for faster resampling
bathymetry_resampled = bathymetry.isel(
    lat=slice(None, None, sampling_factor),
    lon=slice(None, None, sampling_factor)
)

# Plot the data with a colormap
bathymetry_resampled.plot(cmap=cmocean.cm.topo, vmin=vmin, vmax=vmax)

# Plot contours
bathymetry_resampled.plot.contour(levels=contour_levels, colors='black', linewidths=0.5)

plt.show()
