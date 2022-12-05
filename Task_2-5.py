# -*- coding: utf-8 -*-
"""
Author: Ashley Morris
Date: 18/10/2022
Title: Task_2-5


@author:         Ashley Morris
@organization:   University of Tasmania
@detail:         Task 2.5: This program will implement kernel operations on
                 spatial raster layers. Specifically the Fledermaus 'simple'
                 slope calculation. Accepts larger kernel sizes.
                
@created:     18/10/2022
@modified:    18/10/2022
@version:     0.1

@change:     18/10/2022 (0.1) Initial coding
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
output_file_name_a = data_location / "mean_slope_raster_task_2-5.tif"
output_file_name_b = data_location / "min_slope_raster_task_2-5.tif"
output_file_name_c = data_location / "max_slope_raster_task_2-5.tif"

start_time = time.perf_counter()

###########################################################
# GLOBAL CONSTANTS
###########################################################
VERSION = 0.1
OFFSET = 3.86 # N-value for converting to orthometric height
NODATA = -9999 

PLTOUT = True # Set to True if you want to ouput the plots for slope rasters
RASTEROUT = False # Set to True if you want to output the new raster

DEMPROJ = 7855 # I am assuming that the DEM was created using MGA 2020 Zone 55
ORIGIN_X = 50342.2695 # Used as the coordinate for the origin point of the DEM
ORIGIN_Y = 5239398.2560
CELL_SIZE = 0.1 # The DEM cell size in metres

KERNEL_SIZE = 9 # Determines the dimensions of the kernel i.e 5 = 5x5
KERNEL_OFFSET = KERNEL_SIZE//2

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

# Calculate the local means for slope using the Fledermaus cross shaped kernel
def raster_fledermaus_mean_min_max(inarr, nrows, ncols):
    """
    @brief: A nested for-loop with a n x n kernel operator. Determines the mean, min and max slope from the fledermaus kernel.
    @author: Ashley Morris
    @created: 17/10/2022
    @param: inarr: The input data array, nrows: The number of rows, ncols: The number columns
    @return: Returned in the following order - 
        outarr_mean: The output array containng values derived from the mean kernel operation
        outarr_min: The output array containing the minimum slope values from the kernel operation
        outarr_max: The output array containing the maximum slope values from the kernel operation  
    """
    # Creating empty arrays to contain the new output arrays
    outarr_mean = np.zeros((nrows, ncols), dtype = np.float32)
    outarr_min = np.zeros((nrows, ncols), dtype = np.float32)
    outarr_max = np.zeros((nrows, ncols), dtype = np.float32)
    
    #Kernel Operation - using list slicing, making sure to start at an offset that will ensure the full range of the kernel is captured
    for r in range(KERNEL_OFFSET, nrows - KERNEL_OFFSET):
        for c in range(KERNEL_OFFSET, ncols - KERNEL_OFFSET):
            
            # Reading a n x n kernel - the KERNEL_SIZE should be set at the top of script
            kernel = outgrid[r - KERNEL_OFFSET : r + KERNEL_OFFSET + 1, c - KERNEL_OFFSET : c + KERNEL_OFFSET + 1]
            
            # Retrieve the centre value of the cross
            centre_val = kernel[KERNEL_OFFSET : KERNEL_OFFSET + 1, KERNEL_OFFSET : KERNEL_OFFSET + 1]
            c_val = centre_val[0][0]
            
            # Retrieve the western arm of the cross
            west_vals = kernel[KERNEL_OFFSET : KERNEL_OFFSET + 1, : KERNEL_OFFSET]
            
            # Retrieve the eastern arm of the cross
            east_vals = kernel[KERNEL_OFFSET : KERNEL_OFFSET + 1, KERNEL_OFFSET + 1 : ]
            
            # Retrieve the northern arm of the cross
            north_vals = kernel[ : KERNEL_OFFSET, KERNEL_OFFSET : KERNEL_OFFSET + 1]
            
            # Retrieve the southern arm of the cross
            south_vals = kernel[KERNEL_OFFSET + 1 : , KERNEL_OFFSET : KERNEL_OFFSET + 1]
            
            # This list will contain the values of fledermaus cross kernel - this creates a ragged tuple, don't know how to fix this
            #kernel_vals = np.array([north_vals, east_vals, south_vals, west_vals, centre_val], dtype=object)
            
            ##############
            # Calculations
            ##############
            
            # Slope for West Arm
            west_slope = np.arctan((west_vals[0][0] - c_val) / (2 * CELL_SIZE))
            
            # Slope for East Arm
            east_slope = np.arctan((east_vals[0][KERNEL_OFFSET - 1] - c_val) / (2 * CELL_SIZE))
            
            # Slope for North Arm
            north_slope = np.arctan((north_vals[0][0] - c_val) / (2 * CELL_SIZE))

            # Slope of the South Arm
            south_slope = np.arctan((south_vals[KERNEL_OFFSET - 1][0] - c_val) / (2 * CELL_SIZE))            
            
            # The average slope as I interpreted from the spec (*The reported result is the arctan function of the average value of the four slopes, coverted to degrees)
            avg_slope = np.arctan((west_slope + east_slope + north_slope + south_slope) / 4) * (180/np.pi)
            
            # The minimum slope out of the four directions
            min_slope = np.min([west_slope, east_slope, north_slope, south_slope]) * (180/np.pi)
            
            # The maximum slope out of the four directions
            max_slope = np.max([west_slope, east_slope, north_slope, south_slope]) * (180/np.pi)
            
            # Assigning to the output arrays
            outarr_mean[r,c] = avg_slope
            outarr_min[r,c] = min_slope
            outarr_max[r,c] = max_slope
            
    return outarr_mean, outarr_min, outarr_max

def plot_mean_dem(slope_raster):
    """
    @brief: Plots the Slope Mean from the Fledermaus kernel operation.
    @author: Ashley Morris
    @created: 18/10/2022
    @param: slope_raster: Expects the input data array that represents the mean slope raster
    """
    fig = plt.figure(frameon = False)
    im1 = plt.imshow(slope_raster, cmap='Greys') # Image 1 - Slope raster
    cbar = plt.colorbar()
    cbar.set_label('Average Slope (degrees)', rotation = 270, labelpad = 20) # Set colourbar title
    
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('Pixel No. (y)', labelpad = 20)
    plt.xlabel('Pixel No. (x)', labelpad = 20)
    plt.title('Homehill Landslide: Slope Raster (Using Mean Kernel) - Task 2.5')
    plt.show()
    
def plot_min_dem(slope_raster):
    """
    @brief: Plots the Slope Minimum from the Fledermaus kernel operation.
    @author: Ashley Morris
    @created: 18/10/2022
    @param: slope_raster: Expects the input data array that represents the mean slope raster
    """
    fig = plt.figure(frameon = False)
    im1 = plt.imshow(slope_raster, cmap='Greys') # Image 1 - Slope raster
    cbar = plt.colorbar()
    cbar.set_label('Minimum Slope (degrees)', rotation = 270, labelpad = 20) # Set colourbar title
    
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('Pixel No. (y)', labelpad = 20)
    plt.xlabel('Pixel No. (x)', labelpad = 20)
    plt.title('Homehill Landslide: Slope Raster (Using Min Kernel Values) - Task 2.5')
    plt.show()
    
def plot_max_dem(slope_raster):
    """
    @brief: Plots the Slope Maximum from the Fledermaus kernel operation.
    @author: Ashley Morris
    @created: 18/10/2022
    @param: slope_raster: Expects the input data array that represents the mean slope raster
    """
    fig = plt.figure(frameon = False)
    im1 = plt.imshow(slope_raster, cmap='Greys') # Image 1 - Slope raster
    cbar = plt.colorbar()
    cbar.set_label('Maximum Slope (degrees)', rotation = 270, labelpad = 20) # Set colourbar title
    
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('Pixel No. (y)', labelpad = 20)
    plt.xlabel('Pixel No. (x)', labelpad = 20)
    plt.title('Homehill Landslide: Slope Raster (Using Max Kernel Values) - Task 2.5')
    plt.show()

def export_raster_geotiff(raster, output_file_name):
    '''
    @brief: Exports rasters as GeoTiffs.
    @author: Ashley Morris
    @created: 18/10/2022
    @param: raster: The input raster to be converted, output_file_name: The path and file name of the output file
    '''
    driver = "GTiff" # The type of file to create
    dim = raster.shape #dimension of grid
    height = dim[0] # The number of rows
    width = dim[1] # The number of columns
    count = 2      # The number of bands
    dtype = raster.dtype #The data type of the array
    crs = CRS.from_epsg(DEMPROJ) # The coordinate reference system
    transform = from_origin(ORIGIN_X, ORIGIN_Y, CELL_SIZE, CELL_SIZE)
    with rasterio.open(output_file_name, "w",
                       driver = driver,
                       height = height,
                       width = width,
                       count = count,
                       dtype = dtype,
                       crs = crs,
                       transform = transform) as rout:
        rout.write(raster, indexes = 1)
        
##############
# Start Here #
##############

print("\t-----------")
print("\t Task_2-5.py")
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
    
    # Output arrays of the fledermaus kernel operation
    fledermaus_dem_mean, fledermaus_dem_min, fledermaus_dem_max = raster_fledermaus_mean_min_max(outgrid, nrows, ncols)
    
    
    #Plot the 3 different slope layers using matplotlib
    if PLTOUT:
        plot_mean_dem(fledermaus_dem_mean)
        plot_min_dem(fledermaus_dem_min)
        plot_max_dem(fledermaus_dem_max)
    
    # This section will be executed if the RASTEROUT flag is set to True
    # Outputs the slope layers as a GeoTIFF
    if RASTEROUT:
        #Export Mean Slope Raster
        export_raster_geotiff(fledermaus_dem_mean, output_file_name_a)
        
        #Export Min Slope Raster
        export_raster_geotiff(fledermaus_dem_min, output_file_name_b)
        
        #Export Max Slope Raster
        export_raster_geotiff(fledermaus_dem_max, output_file_name_c)

end_time = time.perf_counter()
code_time = end_time - start_time
print("\n")
print("EXECUTION TIME: This script completed in {:.6f} seconds".format(code_time))