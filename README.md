# Blakes-Opening-UAV-Data-Applications
This repository contains several individual python scripts that address particular tasks that make use of GPS data retrieved from UAV flights over Blakes Opening in Tasmania. 

*EXTERNAL LIBRARIES USED:*
pyproj
matplotlib (pyplot, patches)
scipy (spatial, stats)
numpy
rasterio (crs, transform)

*Script 1:*
This program identifies the maximum spatial extent of the UAV flight path (bounding box), as well as find and show the minimum and maximum easting, northing, and height values
of the provided .pos file.

*Script 2:*
This program produces a 2D scatterplot of the UAV flight path, showing the flight path bounding box as a polygon, and convex hull of the flight path as a polygon.

*Script 3:*
This program will read a DEM grid layer and query and manipulate the values. The output will produce a hillshade layer.

*Script 4:*
This program implements kernel operations on spatial raster layers. Specifically the Fledermaus 'simple' slope calculation.

*Script 5:*
This program will implement kernel operations on spatial raster layers. Specifically the Fledermaus 'simple' slope calculation. Accepts larger kernel sizes.

*Script 6:*
This program will implement the 2D Area:3D Area ratio. This ratio aims to quanitfy the roughness of terrain in terms of the planimetric-to-surface area ratio.
This script can take a while to run.

**Date last updated: 05/12/2022**
