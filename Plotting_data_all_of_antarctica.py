#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
from siphon.catalog import TDSCatalog
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import ticker as mticker

# Setting up the figure and related variables

# Set up the map projection (you can choose a different projection if needed)
projection = ccrs.SouthPolarStereo()
transform = ccrs.PlateCarree()

# Create a figure and axis
fig, ax = plt.subplots(subplot_kw={'projection': projection}, figsize=(30, 25))

# Geospatial range to plot
# Full range
lat_min = -90
lat_max = -57
lon_min = -180
lon_max = 180
# Zoom
# lat_min = -80
# lat_max = -70
# lon_min = 0
# lon_max = 60

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
        print(f"Processing dataset {n}: {dataset.name}")
        ds = xr.open_dataset(dataset.access_urls['OPENDAP'])
        bathymetry = ds['z']

        # Select every nth sample for faster resampling
        bathymetry_resampled = bathymetry.isel(
            lat=slice(None, None, sampling_factor),
            lon=slice(None, None, sampling_factor)
        )

        # Selecting data only within geospatial limits specified
        bathymetry_resampled = bathymetry_resampled.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
        if bathymetry_resampled.size == 0:
            print("No data in the specified range for this file.")
            continue  # Skip this file and move to the next one

        # Plot tile
        im =  bathymetry_resampled.plot(cmap=cmocean.cm.topo, vmin=vmin, vmax=vmax, ax=ax, transform=transform, add_colorbar=False)
        # Plot contours
        bathymetry_resampled.plot.contour(levels=contour_levels, colors='black', linewidths=0.1)

# Add labels, title, colorbar, etc. as needed
ax.set_title('Bathymetry of Antarctica below the ice')

# Configure gridlines
gl = ax.gridlines(
    crs=transform, draw_labels=True, linewidth=0.5,
    color='gray', alpha=0.7, linestyle='--'
)
gl.ylocator = mticker.AutoLocator()
gl.xlocator = mticker.AutoLocator()
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}

ax.set_extent([lon_min, lon_max, lat_min, lat_max], transform)

# Create a single colorbar for both plots
cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label('Bathymetry (meters)')

# Show the plot
plt.savefig('antarctica.png', dpi=500)

plt.show()