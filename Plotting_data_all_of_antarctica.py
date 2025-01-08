#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
from siphon.catalog import TDSCatalog

# Setting up the figure and related variables

# Set up the map projection (you can choose a different projection if needed)
projection = ccrs.SouthPolarStereo()

# Create a figure and axis
fig, ax = plt.subplots(subplot_kw={'projection': projection}, figsize=(30, 25))

# Geospatial range to plot
lat_min = -80
lat_max = -70
lon_min = 30
lon_max = 80
# Elevation range for colour scale
vmax = 3000
vmin = vmax * -1
# Contour interval
contour_interval = 500
contour_levels = np.arange(vmin, vmax + contour_interval, contour_interval)
# Plot only every nth sample in both lat and lon to speed up processing
sampling_factor = 10

# Traversing the THREDDS server
catalog_url = 'https://www.ngdc.noaa.gov/thredds/catalog/global/ETOPO2022/15s/15s_bed_elev_netcdf/catalog.xml'

# Access the THREDDS catalog
catalog = TDSCatalog(catalog_url)

# Traverse through the catalog and print a list of the NetCDF files
datasets_filenames = catalog.datasets

# Traverse through the catalog and print URLs of the NetCDF files
datasets_urls = []
for dataset in catalog.datasets.values():
    datasets_urls.append(dataset.access_urls['OPENDAP'])

n=0
for dataset in catalog.datasets.values():
    if 'ETOPO_2022_v1_15s_S' in dataset.name:
        n=n+1
        print(n)
        ds = xr.open_dataset(dataset.access_urls['OPENDAP'])
        bathymetry = ds['z']

        # Selecting data only within geospatial limits specified
        bathymetry_zoomed = bathymetry.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        # Check if there is data in the subset
        if bathymetry_zoomed.size == 0:
            print("No data in the specified range for file")
            continue # Skip this file and move to the next one
        else:
            bathymetry = bathymetry_zoomed

        # Select every nth sample for faster resampling
        bathymetry_resampled = bathymetry.isel(
            lat=slice(None, None, sampling_factor),
            lon=slice(None, None, sampling_factor)
        )
        # Plot tile
        im =  bathymetry_resampled.plot(cmap=cmocean.cm.topo, vmin=vmin, vmax=vmax, ax=ax, transform=ccrs.PlateCarree(), add_colorbar=False)
        # Plot contours
        bathymetry_resampled.plot.contour(levels=contour_levels, colors='black', linewidths=0.1)

# Add labels, title, colorbar, etc. as needed
ax.set_title('Bathymetry Comparison')
#ax.set_xlabel('Longitude')
#ax.set_ylabel('Latitude')

ax.set_extent([-180, 180, -90, -57], ccrs.PlateCarree())

# Create a single colorbar for both plots
cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label('Bathymetry (meters)')

# Show the plot
plt.savefig('antarctica.png', dpi=500)

plt.show()