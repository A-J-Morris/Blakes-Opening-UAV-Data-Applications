# -*- coding: utf-8 -*-
"""
Author: Ashley Morris
Date: 17/10/2022
Title: Task_2-3


@author:         Ashley Morris
@organization:   University of Tasmania
@detail:         Task 2.3: This program will read a DEM grid layer
                 and query and manipulate the values. The output will produce
                 a hillshade layer.
                
@created:     17/10/2022
@modified:    17/10/2022
@version:     0.1

@change:     17/10/2022 (0.1) Initial coding
"""

###########################################################
# PACKAGE & LIBRARY IMPORTS
###########################################################

# Inbuilt Libraries:
from pathlib import Path
import time
import warnings

#External Libraries:
import matplotlib.pyplot as plt
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
import numpy as np

###########################################################
# USER DEFINED VARIABLES
###########################################################

#Set this to the location where the data is stored on YOUR device
data_location = Path("F:/UTAS_2022/Advanced GIS/Assignment 3 Programming")

#Make sure the NAME of the FILE you are trying to access is CORRECT 
file_name = data_location / "Homehill_2020_DSM_30cm.tif"
file = (data_location / file_name).as_posix()
output_file_name = data_location / "outputraster.tif"

start_time = time.perf_counter()

###########################################################
# GLOBAL CONSTANTS
###########################################################
VERSION = 0.1
OFFSET = 3.86 # N-value for converting to orthometric height
NODATA = -9999 

RASTEROUT = False # Set to True if you want to output the new raster
DEMPROJ = 7855 # I am assuming that the DEM was created using MGA 2020 Zone 55
ORIGIN_X = 50342.2695 # Used as the coordinate for the origin point of the DEM
ORIGIN_Y = 5239398.2560
CELL_SIZE = 0.1 # The DEM cell size in metres

# Azimuth and elevation in degrees for calculating hillshade
AZIMUTH_DEG = 340
ZENITH_DEG = 60

###########################################################
# USER DEFINED FUNCTIONS
###########################################################

def raster_add(inarr, OFFSET, nrows, ncols):
    """
    @brief: Add a value to a grid
    @author: Ashley Morris
    @created: 17/10/2022

    Parameters
    ----------
    inarr : LIST
        The input data array.
    OFFSET : FLOAT
        The height offset value for converting to orthometric height.
    nrows : INT
        The number of rows.
    ncols : INT
        The number of columns.

    Returns
    -------
    @return: newarr: The new output array based on the modified input array.
    """
    
    outarr = inarr + OFFSET
    
    return outarr

def hillshade (arr, azimuth, angle_altitude):
    '''
    @brief: Generate a hillshade layer from a given DEM (as an array)
    @source: This function was developed by Roger Veciana, and can be found\
              at https://github.com/rveciana/introduccion-python-geoespacial/blob/master/hillshade.py
    @params: arr: The DEM as an array, azimuth: The azimuth of the illumination source, angle_altitude: The zenith (altitude)
    @return: A hypothetical illumination value for each pixel based on provided parameters
    '''
    azimuth = 360.0 - azimuth
    
    x, y = np.gradient(arr)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi/180.
    altituderad = angle_altitude*np.pi/180.
    
    shaded = np.sin(altituderad)*np.sin(slope) + np.cos(altituderad)*np.cos(slope)*np.cos((azimuthrad - np.pi/2.) - aspect)
    
    return 255*(shaded + 1)/2
    
##############
# Start Here #
##############

print("\t-----------")
print("\t Task_2-3.py")
print("\t-----------")
print("\t Version " + str(VERSION))

###########################################################
# MAIN SCRIPT
###########################################################
warnings.filterwarnings('ignore')

# Reading a DEM grid layer as a GeoTIFF file using rasterio
with rasterio.open(file) as dataset:
    dem = dataset.read(1, masked = True) # Reading as a masked array to deal with nodata values
    nrows, ncols = dem.shape
    outgrid = raster_add(dem, OFFSET, nrows, ncols)
    
    # Replace any numpy nodata values with -9999 for output (just in case it occurs)
    outgrid[np.isnan(outgrid)] = NODATA
    
    start_time_hillshade = time.perf_counter() # Start the timer for hillshade process
    
    # Generate the hillshade array
    hs_array = hillshade(outgrid, AZIMUTH_DEG, ZENITH_DEG)
    
    end_time_hillshade = time.perf_counter() # End the timer for hillshade process
    
    # This section will be executed if the RASTEROUT flag is set to True
    # Outputs the hillshade raster layer as a GeoTIFF
    if RASTEROUT:
        driver = "GTiff" # The type of file to create
        dim = hs_array.shape
        height = dim[0] # The number of rows
        width = dim[1] # The number of columns
        count = 2      # The number of bands
        dtype = hs_array.dtype #The data type of the array
        crs = CRS.from_epsg(DEMPROJ) # The coordinate reference system
        transform = from_origin(ORIGIN_X, ORIGIN_Y, CELL_SIZE, CELL_SIZE)
        with rasterio.open(output_file_name, "w",
                           driver = driver,
                           height = height,
                           width = width,
                           count = count,
                           dtype = dtype,
                           crs = crs,
                           transform = transform) as hshade:
            hshade.write(hs_array, indexes = 1)
    
    fig = plt.figure(frameon = False)
    im1 = plt.imshow(outgrid, cmap='terrain_r') # Image 1 - DEM
    cbar = plt.colorbar()
    cbar.set_label('Elevation (m)', rotation = 270, labelpad = 20) # Set colourbar title
    
    im2 = plt.imshow(hs_array, cmap = 'Greys', alpha = 0.8) # Image 2 - Hillshade
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('Pixel No. (y)', labelpad = 20)
    plt.xlabel('Pixel No. (x)', labelpad = 20)
    plt.title('HomeHill Landslide: DEM + Hillshade')
    plt.show()

end_time = time.perf_counter()
code_time = end_time - start_time
hillshade_time = end_time_hillshade - start_time_hillshade
print("\n")
print("EXECUTION TIME: This script completed in {:.6f} seconds".format(code_time))
print("HILLSHADE TIME: The hillshade process completed in {:.6f} seconds".format(hillshade_time))