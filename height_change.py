# height_change.py
# by Tyler Wixtrom and Kevin Goebbert
# This script will change WRF geogrid points defined to a new terrain height.
# It will output a new version of the geo_em.dxx.nc file created by WPS
# This can be used to replace the file found in the WRF/WPS directory.
# Repeating WRF startup steps after execution of geogrid.exe will allow the
# model to be run with the desired water region defined as land.
#
#
# The region of the geography grid that you wish to modify should be defined by
# the lat/lon coordinates of the corners. Any gridpoints within this region with
# the landmask set to water (landmask = 0) will be modified with the attributes
# of the point defined in part b
#
# Important information about what fields needed to be modified came from the
# following website: http://www.wrfems.info/viewtopic.php?t=141
# Code for copying netCDF variables and dimensions from one file to another in
# Python came from the following website:
# http://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one


# import libraries
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt

from netCDF4 import Dataset

import numpy as np
import numpy.ma as ma

# Part A: Define location of grid
# Define lower left and upper right corner of area to change to land
UR = [48.7, -83.7]
UL = [48.7, -92.3]
LL = [46.65, -92.3]
LR = [46.65, -83.7]


# Part B: set point for new data
# lat/lon of location to use for the new land properties over the lake
# the lake will be redefined to match the gridpoint nearest to the point defined below
llat = 47.56
llon = -91.30

# This will make a copy of the input file for editing
# The file will be written after data is edited
# code from Kevin Goebbert

# input file
print('opening input file...')
dsin = Dataset('in.nc', 'r')

# output file
dsout = Dataset('out.nc', 'w')

# Copy dimensions
for dname, the_dim in dsin.dimensions.items():
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

# Pulling just the latitude and londitude values from the computational mass field
# lat and lon grids
# NOTE: Lats and Lons are two dimensional as stored in geo_em.d01.nc
# lat.shape -> (1,104,124) for grid I am using (Valpo WRF run)
# This will make finding a subset area a little harder...
latA = dsin.variables['CLAT']
lonA = dsin.variables['CLONG']
lat = latA[0, :, :]  # Reducing dimensionality by 1 (getting rid of time dimension)
lon = lonA[0, :, :]

# Find lat_lon_index for a 2D lat/lon array


def lat_lon_2D_index(y, x, lat1, lon1):
    """
    This function calculates the distance from a desired lat/lon point
    to each element of a 2D array of lat/lon values, typically from model output,
    and determines the index value corresponding to the nearest lat/lon grid point.

    x = longitude array
    y = latitude array
    lon1 = longitude point (signle value)
    lat1 = latitude point (single value)

    Returns the index value for nearest lat/lon point on grid

    Equations for variable distiance between longitudes from
    http://andrew.hedges.name/experiments/haversine/
    """
    R = 6373.*1000.  # Earth's Radius in meters
    rad = np.pi/180.
    x1 = np.ones(x.shape)*lon1
    y1 = np.ones(y.shape)*lat1
    dlon = np.abs(x-x1)
    dlat = np.abs(y-y1)
    a = (np.sin(rad*dlat/2.))**2 + np.cos(rad*y1) * np.cos(rad*y) * (np.sin(rad*dlon/2.))**2
    c = 2 * np.arctan2( np.sqrt(a), np.sqrt(1-a))
    d = R * c
    return np.unravel_index(d.argmin(), d.shape)


# Finding the index values for the UR and LL box to limit the area where the landmask will be changed
ilat_UR, ilon_UR = lat_lon_2D_index(lat, lon, UR[0], UR[1])
ilat_LL, ilon_LL = lat_lon_2D_index(lat, lon, LL[0], LL[1])

# Defining variables to change

# hght_m variable
hgt_m_in = dsin.variables['HGT_M'][:]

mask = ma.make_mask_none(lat.shape)
for i in range(ilat_LL, ilat_UR+1):
    for j in range(ilon_LL, ilon_UR+1):
        mask[i, j] = True

# Getting the index values of gridpoint to be cloned to lake
ilat_lnd, ilon_lnd = lat_lon_2D_index(lat, lon, llat, llon)

# Changing the variables based on above lat/lon point defined
hgt_m_out = hgt_m_in.copy()
hgt_m_out[0, mask == True] = 1500.0

print('writing to output file...')

# Copy variables
for v_name, varin in dsin.variables.items():
    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)

    # Copy variable attributes
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

    if v_name == 'HGT_M':
        outVar[:] = hgt_m_out[:]
    else:
        outVar[:] = varin[:]

# Copy global attributes from infile to outfile
dsout.setncatts({k: dsin.getncattr(k) for k in dsin.ncattrs()})

# Close outfile to write global attributes
dsout.close()

print('Completed Successfully')
