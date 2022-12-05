# -*- coding: utf-8 -*-
"""
Author: Ashley Morris
Date: 17/10/2022
Title: Task_2-2


@author:         Ashley Morris
@organization:   University of Tasmania
@detail:         Task 2.2: This program produces a 2D scatterplot of the UAV
                 flight path, showing the fight path bounding box as a polygon,
                 and convex hull of the flight path as a polygon.
                
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

#External Libraries:
import pyproj
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import ConvexHull
import numpy as np

###########################################################
# USER DEFINED VARIABLES
###########################################################

#Set this to the location where the data is stored on YOUR device
data_location = Path("F:/UTAS_2022/Advanced GIS/Assignment 3 Programming")

#Make sure the NAME of the FILE you are trying to access is CORRECT 
file_name = data_location / '20200910_BlakesOpening_site3_flight1z.pos'
csv_file_name = data_location / 'BlakesOpening_Site3_Flight1z_GDA2020_55.csv'

start_time = time.perf_counter()
easting_list = []
northing_list = []
height_list = []
quality_factor_list = []

###########################################################
# GLOBAL CONSTANTS
###########################################################

#EPSG coordinate systems defined here:
#source of codes available from: https://epsg.io/
INPROJ = pyproj.CRS("epsg:4326") # WGS 84
OUTPROJ = pyproj.CRS("epsg:7855") # GDA2020 MGA Zone 55
CSVOUT = False # Set to True if you wish to output a CSV file
PLTOUT_A = False  # Set to True if you wish to output a scatterplot (A) using Matplotlib
PLTOUT_B = False # Set to True if you wish to output a scatterplot (B) using Matplotlib
PLTOUT_C = False # Set to True if you wish to ouput a scatterplot (C) using Matplotlib
VERSION = 0.1
OFFSET = 3.86 # N-value for converting to orthometric height

#This script will skip this number of HEADER LINES from the input file. Set this to the line at which the tabulated data is provided. 
HEADER_LINES_SKIPPED = 38

#This script will read every STEPth line of the input file. Set to 1 to read ALL lines.
STEP = 10

###########################################################
# USER DEFINED FUNCTIONS
###########################################################

def read_pos_coords(pos_file_line):
    """
    @brief:     Extract coordinate and quality data from a single line within .pos file
    @note:      Outputs in the order X, Y, Z, Q
    @author:    Ashley Morris
    @created:   16/10/2022
    
    Parameters
    ----------
    pos_file_line : STRING
        POS formatted line from the .pos file.

    Returns
    -------
    @return:    lon_dec, lat_dec, height, quality:     Longitude, Latitude, Height (in WGS84), Quality Factor (Q)
    """
    
    line_list = pos_file_line.split(',')
    lat_dec = float(line_list[1]) # Retrieving the signed latitude coordinate in decimal degrees
    lon_dec = float(line_list[2]) # Retrieving the signed longitude coordinate in decimal degrees
    
    # Retrieving the height above the current ellipsoid (metres) -
    # and converting to orthometric height by subtracting the interpolated geoid height.
    height = float(line_list[3]) - OFFSET 
    quality = int(line_list[4]) # Retrieving the quality factor (1-6)
    
    return lon_dec, lat_dec, height, quality 

def calc_bounding_box(x_list, y_list):
    """
    @brief:     Calculate the bounding box
    @note:      Works with cartesian coordinates and longitude / Latitude
    @author:    Ashley Morris
    @created:   16/10/2022
    
    Parameters
    ----------
    x_list : LIST
        Expecting a list of easting coordinates.

    y_list : LIST
        Expecting a list of northing coordinates.
        
    Returns
    -------
    @return:    x_min, x_max, y_min, y_max, width, height:     Outputs the bounding box coordinates, width, and height.
    """
    
    x_min = min(x_list)
    y_min = min(y_list)
    x_max = max(x_list)
    y_max = max(y_list)
    
    width = x_max - x_min
    height = y_max - y_min
    
    return x_min, x_max, y_min, y_max, width, height

##############
# Start Here #
##############

print("\t-----------")
print("\t Task_2-2.py")
print("\t-----------")
print("\t Version " + str(VERSION))

###########################################################
# MAIN SCRIPT
###########################################################
with open(file_name, 'r') as f:
    lines = f.readlines()
    num_lines = len(lines)
    
    for line_num in range(HEADER_LINES_SKIPPED, num_lines, STEP):
        line = lines[line_num]
        
        #Extract the coordinate, height, and quality data
        lon_dec, lat_dec, height, quality = read_pos_coords(line)
        
        #Project to GDA2020 MGA Zone 55
        transformer = pyproj.Transformer.from_crs(INPROJ, OUTPROJ, always_xy=True)
        
        #Always check the correct GIS input order for X and Y:
        easting, northing = transformer.transform(lon_dec, lat_dec)
        
        easting_list.append(easting)
        northing_list.append(northing)
        height_list.append(height)
        quality_factor_list.append(quality)

# Retrieve the bounding box attributes
x_min, x_max, y_min, y_max, width, height = calc_bounding_box(easting_list, northing_list)

###########################################################
# (OPTIONAL - SET FLAG) CREATE PLOTS
###########################################################
if PLTOUT_A:
    x = easting_list
    y = northing_list
    
    fig, ax = plt.subplots()
    
    sc = ax.scatter(x, y, c = height_list, alpha = 0.5) # Points are coloured by height
    
    # Setting the axis and title labels
    plt.title('Blakes Opening: UAV Flight Path')
    ax.set_ylabel("Northings (m)")
    ax.set_xlabel("Eastings (m)")
    
    cbar = fig.colorbar(sc)
    cbar.set_label("Orthometric Height (m)", loc='center')
    plt.ticklabel_format(style='plain')
    plt.show()

if PLTOUT_B:
    x = easting_list
    y = northing_list
    
    fig, ax = plt.subplots()
    
    sc = ax.scatter(x, y, c = height_list, alpha = 0.5) # Points are coloured by height
    
    # Adding the bounding box for the flight path as a rectangular polygon using bounding box calculations
    ax.add_patch (
        patches.Rectangle(
            xy = (x_min, y_min),
            width = width,
            height = height,
            linewidth = 1,
            color = 'red',
            fill = False
        )
    )
    
    # Setting the axis and title labels
    plt.title('Blakes Opening: UAV Flight Path')
    ax.set_ylabel("Northings (m)")
    ax.set_xlabel("Eastings (m)")
    
    cbar = fig.colorbar(sc)
    cbar.set_label("Orthometric Height (m)", loc='center')
    plt.ticklabel_format(style='plain')
    plt.show()
    
if PLTOUT_C:
    x = easting_list
    y = northing_list
    
    fig, ax = plt.subplots()
    
    sc = ax.scatter(x, y, c = height_list, alpha = 0.5) # Points are coloured by height
    
    points = np.stack((x, y), axis = 1) # Combine the eastings and northings into one list of coordinate points
    hull = ConvexHull(points)
    
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k--')
    
    # Setting the axis and title labels
    plt.title('Blakes Opening: UAV Flight Path (With Convex Hull)')
    ax.set_ylabel("Northings (m)")
    ax.set_xlabel("Eastings (m)")
    
    cbar = fig.colorbar(sc)
    cbar.set_label("Orthometric Height (m)", loc='center')
    plt.ticklabel_format(style='plain')
    plt.show()
    
###########################################################
# (OPTIONAL - SET FLAG) WRITE CSV FILE
###########################################################
if CSVOUT:
    with open(csv_file_name, 'w') as csvout:
        csvout.write("Point ID, Easting, Northing, Height, Quality\n")
        num_coords = len(easting_list)

        for i in range(num_coords):
            point_id = (i * STEP) + 1
                
            csvout.write("{},{:.3f},{:.3f},{:.3f},{}\n".format(point_id, easting_list[i], northing_list[i], height_list[i], quality_factor_list[i]))

end_time = time.perf_counter()
code_time = end_time - start_time
print("\n")
print("EXECUTION TIME: This script completed in {:.6f} seconds".format(code_time))