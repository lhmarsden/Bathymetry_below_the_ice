#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
from siphon.catalog import TDSCatalog
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import ticker as mticker

# Set up the map projection
projection = ccrs.NorthPolarStereo(central_longitude=-45)
transform = ccrs.PlateCarree()

# Create a figure and axis
fig, ax = plt.subplots(subplot_kw={'projection': projection}, figsize=(12, 10))

# Geospatial range to plot
# Set to 'False' to plot full range of the data, or provide a value
zoom = True
if zoom is True:
    lat_min = 58
    lat_max = 85
    lon_min = -80
    lon_max = -10
else:
    lat_min = None
    lat_max = None
    lon_min = None
    lon_max = None

# Initialising values
computed_lat_min = float('inf')
computed_lat_max = float('-inf')
computed_lon_min = float('inf')
computed_lon_max = float('-inf')

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

n = 0
for dataset in catalog.datasets.values():
    if 'ETOPO_2022_v1_15s_N' in dataset.name:
        n += 1
        print(f"Processing dataset {n}: {dataset.name}")
        ds = xr.open_dataset(dataset.access_urls['OPENDAP'])
        bathymetry = ds['z']
        bathymetry_resampled = bathymetry.isel(
            lat=slice(None, None, sampling_factor),
            lon=slice(None, None, sampling_factor)
        )

        if zoom == True:
            # Selecting data only within geospatial limits specified
            bathymetry_resampled = bathymetry_resampled.sel(lat=slice(lat_min, lat_max), lon=slice(lon_min, lon_max))
            if bathymetry_resampled.size == 0:
                print("No data in the specified range for this file.")
                continue  # Skip this file and move to the next one

        # Update the global lat_min, lat_max, lon_min, lon_max across all files
        if zoom is False:
            computed_lat_min = min(computed_lat_min, bathymetry_resampled.coords['lat'].min().values)
            computed_lat_max = max(computed_lat_max, bathymetry_resampled.coords['lat'].max().values)
            computed_lon_min = min(computed_lon_min, bathymetry_resampled.coords['lon'].min().values)
            computed_lon_max = max(computed_lon_max, bathymetry_resampled.coords['lon'].max().values)

        # Plot the data
        im = bathymetry_resampled.plot(
            cmap=cmocean.cm.topo, vmin=vmin, vmax=vmax,
            ax=ax, transform=transform, add_colorbar=False
        )
        # Plot contours
        bathymetry_resampled.plot.contour(
            ax=ax, levels=contour_levels, colors='black',
            linewidths=0.5, transform=transform
        )

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

# Clip the map to the data extent
if zoom is False:
    lat_min = computed_lat_min
    lat_max = computed_lat_max
    lon_min = computed_lon_min
    lon_max = computed_lon_max
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=transform)

# Add title and colorbar
ax.set_title('Bathymetry of Greenland below the ice', fontsize=14)
cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.1)
cbar.set_label('Bathymetry (meters)', fontsize=12)

# Save the plot
plt.savefig('greenland.png', dpi=500)

# Show the plot
plt.show()