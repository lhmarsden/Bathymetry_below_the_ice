#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean

# Traverse THREDDS server, isolate 'N' for Greenland.
url = 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/15s/15s_bed_elev_netcdf/ETOPO_2022_v1_15s_N75W030_bed.nc'

ds = xr.open_dataset(url)

bathymetry = ds['z']

# Select every nth sample for faster resampling
sampling_factor = 10
bathymetry_resampled = bathymetry.isel(
    lat=slice(None, None, sampling_factor),
    lon=slice(None, None, sampling_factor)
)

# Define plot parameters
vmax = 3000
vmin = vmax * -1

# Plot the data with a colormap
bathymetry_resampled.plot(cmap=cmocean.cm.topo, vmin=vmin, vmax=vmax)

# Plot contours
contour_interval = 500
contour_levels = np.arange(vmin, vmax + contour_interval, contour_interval)
bathymetry_resampled.plot.contour(levels=contour_levels, colors='black', linewidths=0.5)

plt.show()
