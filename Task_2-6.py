# -*- coding: utf-8 -*-
"""
Author: Ashley Morris
Date: 19/10/2022
Title: Task_2-6


@author:         Ashley Morris
@organization:   University of Tasmania
@detail:         Task 2.6: This program will implement the 2D Area:3D Area Ratio.
                           This ratio aims to quantify the roughness of terrain in
                           terms of the planimetric-to-surface area ratio. The script can take a while to run.
@created:     19/10/2022
@modified:    19/10/2022
@version:     0.1

@change:     19/10/2022 (0.1) Initial coding
"""

###########################################################
# PACKAGE & LIBRARY IMPORTS
###########################################################

# Inbuilt Libraries:
from pathlib import Path
import time
import warnings
import math

#External Libraries:
import matplotlib.pyplot as plt
import rasterio
from rasterio.crs import CRS
from rasterio.transform import from_origin
from scipy.stats import kurtosis
import numpy as np

###########################################################
# USER DEFINED VARIABLES
###########################################################

#Set this to the location where the data is stored on YOUR device
data_location = Path("F:/UTAS_2022/Advanced GIS/Assignment 3 Programming")

#Make sure the NAME of the FILE you are trying to access is CORRECT 
file_name = data_location / "Homehill_2020_DSM_30cm.tif"
file = (data_location / file_name).as_posix()
output_file_name_a = data_location / "raster_task_2-6.tif"

start_time = time.perf_counter()

###########################################################
# GLOBAL CONSTANTS
###########################################################
VERSION = 0.1
OFFSET = 3.86 # N-value for converting to orthometric height
NODATA = -9999 

PLTOUT = False # Set to True if you want to ouput the plots for slope rasters
RASTEROUT = False # Set to True if you want to output the new raster

DEMPROJ = 7855 # I am assuming that the DEM was created using MGA 2020 Zone 55
ORIGIN_X = 50342.2695 # Used as the coordinate for the origin point of the DEM
ORIGIN_Y = 5239398.2560
CELL_SIZE = 0.1 # The DEM cell size in metres

KERNEL_SIZE = 3 # Determines the dimensions of the kernel i.e 5 = 5x5 (3x3 is used for area ratio index)
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

def raster_area_ratio_index(inarr, nrows, ncols):
    """
    @brief: Uses the 2D Area: 3D Area ratio index calculation from Jenness (2004) to calculate the surface area for each cell. NOTE: METHOD IS SLOW
    @author: Ashley Morris
    @created: 19/10/2022
    @params: inarr: The input raster as an array, nrows: The number of rows in the input raster, ncols: The number of columns in the input raster
    @return: outarr: The ouput array containing the calculated surface area for each cell from the input grid
    """
    # Creating empty array to contain the new output array
    outarr = np.zeros((nrows, ncols), dtype = np.float32)

    #Kernel Operation - using list slicing, making sure to start at an offset that will ensure the full range of the kernel is captured
    for r in range(KERNEL_OFFSET, nrows - KERNEL_OFFSET):
        for c in range(KERNEL_OFFSET, ncols - KERNEL_OFFSET):
            
            # Reading a n x n kernel - the KERNEL_SIZE should be set at the top of script
            kernel = outgrid[r - KERNEL_OFFSET : r + KERNEL_OFFSET + 1, c - KERNEL_OFFSET : c + KERNEL_OFFSET + 1]
            
            #####################################################################################
            # Area Ratio Index Calculations (using variable names as suggested by Jenness (2004))
            #####################################################################################
            # I have decided to explicitly lay out the steps for ease of following even though they may be able to be completed in less lines
            # and less calculations. It is easier to follow the steps by doing it this way. The assignment DOES NOT ASK for a fast implementation. 
            
            #Getting the value of 'b' (difference in elevation) for each of the two cell pairs (should be 16 edges for 3*3 window)
            # North = n, North-west = nw, North-east = ne, West = w, South-west = sw, South = s, South-east = se, East = e, Central Cell = cent
            
            # Getting the cell values
            centre = kernel[1][1]
            nw = kernel[0][0]
            n = kernel[0][1]
            w = kernel[1][0]
            ne = kernel[0][2]
            sw = kernel[2][0]
            s = kernel[2][1]
            se = kernel[2][2]
            e = kernel[1][2]
            
            # 'b' is the absoulte difference in elevation between the two cells
            b_nw_n = abs(nw - n)
            b_nw_w = abs(nw - w)
            b_nw_cent = abs(nw - centre)
            b_ne_cent = abs(ne - centre)
            b_n_ne = abs(n - ne)
            b_n_cent = abs(n - centre)
            b_w_cent = abs(w - centre)
            b_w_sw = abs(w - sw)
            b_sw_cent = abs(sw - centre)
            b_sw_s = abs(sw - s)
            b_s_cent = abs(s - centre)
            b_s_se = abs(s - se)
            b_se_cent = abs(se - centre)
            b_se_e = abs(se - e)
            b_e_cent = abs(e - centre)
            b_e_ne = abs(e - ne)
            
            # Getting the value for 'a' (planimetric distance) for each cell pair. It is the length of the side of the cells.
            a_e_cent = CELL_SIZE
            a_s_cent = CELL_SIZE
            a_w_cent = CELL_SIZE
            a_n_cent = CELL_SIZE
            a_nw_n = CELL_SIZE
            a_nw_w = CELL_SIZE
            a_w_sw = CELL_SIZE
            a_n_ne = CELL_SIZE
            a_sw_s = CELL_SIZE
            a_s_se = CELL_SIZE
            a_se_e = CELL_SIZE
            a_e_ne = CELL_SIZE
            
            #The planimetric distance must be calculated differently for diagonal cell pairs
            a_ne_cent = math.sqrt(2*(CELL_SIZE ** 2))
            a_se_cent = math.sqrt(2*(CELL_SIZE ** 2))
            a_nw_cent = math.sqrt(2*(CELL_SIZE ** 2))
            a_sw_cent = math.sqrt(2*(CELL_SIZE ** 2))
            
            # Getting the value for 'c' (surface distance) for each cell pair using pythagoras theorem
            # These distances represent the edges of the 8 triangles, but need to be divided by '2' to account for triangles extending past the cell boundaries.
            c_nw_n = math.sqrt((a_nw_n ** 2) + (b_nw_n ** 2)) / 2
            c_nw_w = math.sqrt((a_nw_w ** 2) + (b_nw_w ** 2)) / 2
            c_nw_cent = math.sqrt((a_nw_cent ** 2) + (b_nw_cent ** 2)) / 2
            c_ne_cent = math.sqrt((a_ne_cent ** 2) + (b_ne_cent ** 2)) / 2
            c_n_ne = math.sqrt((a_n_ne ** 2) + (b_n_ne ** 2)) / 2
            c_n_cent = math.sqrt((a_n_cent ** 2) + (b_n_cent ** 2)) / 2
            c_w_cent = math.sqrt((a_w_cent ** 2) + (b_w_cent ** 2)) / 2
            c_w_sw = math.sqrt((a_w_sw ** 2) + (b_w_sw ** 2)) / 2
            c_sw_cent = math.sqrt((a_sw_cent ** 2) + (b_sw_cent ** 2)) / 2
            c_sw_s = math.sqrt((a_sw_s ** 2) + (b_sw_s ** 2)) / 2
            c_s_cent = math.sqrt((a_s_cent ** 2) + (b_s_cent ** 2)) / 2
            c_s_se = math.sqrt((a_s_se ** 2) + (b_s_se ** 2)) / 2
            c_se_cent = math.sqrt((a_se_cent ** 2) + (b_se_cent ** 2)) / 2
            c_se_e = math.sqrt((a_se_e ** 2) + (b_se_e ** 2)) / 2
            c_e_cent = math.sqrt((a_e_cent ** 2) + (b_e_cent ** 2)) / 2
            c_e_ne = math.sqrt((a_e_ne ** 2) + (b_e_ne ** 2)) / 2
            
            # Calculating the areas of the 8 triangles
            
            #Triangle 1 (edges): NW_W, NW_Centre, W_Centre
            s_1 = (c_nw_w + c_nw_cent + c_w_cent) / 2
            tri_1_area = math.sqrt(s_1 * (s_1 - c_nw_w) * (s_1 - c_nw_cent) * (s_1 - c_w_cent))
            
            #Triangle 2 (edges): NW_N, NW_Centre, N_Centre
            s_2 = (c_nw_n + c_nw_cent + c_n_cent) / 2
            tri_2_area = math.sqrt(s_2 * (s_2 - c_nw_n) * (s_2 - c_nw_cent) * (s_2 - c_n_cent))
            
            #Triangle 3 (edges): N_Centre, N_NE, NE_Centre
            s_3 = (c_n_ne + c_ne_cent + c_n_cent) / 2
            tri_3_area = math.sqrt(s_3 * (s_3 - c_n_ne) * (s_3 - c_ne_cent) * (s_3 - c_n_cent))
            
            #Triangle 4 (edges): NE_Centre, E_NE, E_Centre
            s_4 = (c_ne_cent + c_e_ne + c_e_cent) / 2
            tri_4_area = math.sqrt(s_4 * (s_4 - c_ne_cent) * (s_4 - c_e_ne) * (s_4 - c_e_cent))
            
            #Triangle 5 (edges): E_Centre, SE_E, SE_Centre
            s_5 = (c_e_cent + c_se_e + c_se_cent) / 2
            tri_5_area = math.sqrt(s_5 * (s_5 - c_e_cent) * (s_5 - c_se_e) * (s_5 - c_se_cent))
            
            #Triangle 6 (edges): SE_Centre, S_Centre, S_SE
            s_6 = (c_se_cent + c_s_cent + c_s_se) / 2
            tri_6_area = math.sqrt(s_6 * (s_6 - c_se_cent) * (s_6 - c_s_cent) * (s_6 - c_s_se))
            
            #Triangle 7 (edges): S_Centre, SW_Centre, SW_S
            s_7 = (c_s_cent + c_sw_cent + c_sw_s) / 2
            tri_7_area = math.sqrt(s_7 * (s_7 - c_s_cent) * (s_7 - c_sw_cent) * (s_7 - c_sw_s))
            
            #Triangle 8 (edges): W_Centre, SW_Centre, W_SW
            s_8 = (c_w_cent + c_sw_cent + c_w_sw) / 2
            tri_8_area = math.sqrt(s_8 * (s_8 - c_w_cent) * (s_8 - c_sw_cent) * (s_8 - c_w_sw))
            
            # Calculate the TOTAL SURFACE AREA
            total_surface_area = tri_1_area + tri_2_area + tri_3_area + tri_4_area + tri_5_area + tri_6_area + tri_7_area + tri_8_area
            
            # Assigning to the output arrays
            outarr[r,c] = total_surface_area
    
    return outarr

def plot_area_ratio_index(area_ratio_raster):
    """
    @brief: Plots the 2D Area:3D Area ratio layer.
    @author: Ashley Morris
    @created: 19/10/2022
    @param: area_ratio_raster: Expects the input data array that represents 2D Area:3D Area ratio raster layer
    """
    
    fig = plt.figure(frameon = False)
    im1 = plt.imshow(area_ratio_raster, cmap='viridis') # Image 1 - Slope raster
    cbar = plt.colorbar()
    cbar.set_label('Surface Area Derived From 2D:3D Area Ratio (m2)', rotation = 270, labelpad = 20) # Set colourbar title
    
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('Pixel No. (y)', labelpad = 20)
    plt.xlabel('Pixel No. (x)', labelpad = 20)
    plt.title('Homehill Landslide: 2D Area: 3D Area Ratio Index - Task 2.6')
    plt.show()
    
def plot_area_ratio_index_histogram(area_ratio_raster):
    """
    @brief: Plots the 2D Area:3D Area ratio layer as a histogram.
    @author: Ashley Morris
    @created: 19/10/2022
    @param: area_ratio_raster: Expects the input data array that represents 2D Area:3D Area ratio raster layer
    """
    
    fig = plt.figure(frameon = False)
    
    plt.hist(area_ratio_raster, bins = 15)
    
    ax = plt.gca()
    ax.ticklabel_format(useOffset = False, style = 'plain') # Do not use scientific notation
    rotate_x_labels = plt.setp(ax.get_xticklabels(), rotation = 90) # Rotate the x tick labels by 90 degrees
    
    plt.ylabel('No. Observations', labelpad = 20)
    plt.xlabel('Surface Area (m2)', labelpad = 20)
    plt.title('Homehill Landslide: 2D Area: 3D Area Ratio Index - Histogram - Task 2.6')
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
print("\t Task_2-6.py")
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
    
    # Generate the surface area raster using the 2D Area:3D Area ratio index calculations
    surface_area_raster = raster_area_ratio_index(outgrid, nrows, ncols)
    
    #Plot the surface area raster and histogram layers using matplotlib
    if PLTOUT:
        plot_area_ratio_index(surface_area_raster)
        plot_area_ratio_index_histogram(surface_area_raster)
    
    # This section will be executed if the RASTEROUT flag is set to True
    # Outputs the slope layers as a GeoTIFF
    if RASTEROUT:
        #Export 2D:3D Area Ratio layer
        export_raster_geotiff(surface_area_raster, output_file_name_a)
        
        
    # Print descriptive statistics for the 2D Area:3D Area ratio index layer (this is what I interpreted the question to ask for, not the DEM)
    # nanmin, nanmax should exclude and NaN values 
    min_surf_area = np.nanmin(surface_area_raster, axis = None)
    max_surf_area = np.nanmax(surface_area_raster, axis = None)
    mean_surf_area = float(np.nanmean(surface_area_raster))
    stdev_surf_area = float(np.nanstd(surface_area_raster))
    ksis = kurtosis(surface_area_raster, axis = None, nan_policy='omit')
    
    print("\nMinimum Surface Area: {:.6f}\nMaximum Surface Area: {:.6f}\nMean Surface Area: {:.5f}\nSTDEV Surface Area: {:.5f}\nKurtosis: {:.5f}".format(min_surf_area, max_surf_area, mean_surf_area, stdev_surf_area, ksis))
        
end_time = time.perf_counter()
code_time = end_time - start_time
print("\n")
print("EXECUTION TIME: This script completed in {:.6f} seconds".format(code_time))