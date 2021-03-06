# modify_WRF_geogrid

Code for the modification of the WRF-ARW geogrid using Python

by Tyler Wixtrom and Dr. Kevin Goebbert

This code will allow for the removal of water bodies from the WRF geography grid.

To run the WRF without a water body, complete the following steps:
  1. Define region of interest in code. This is the area containing the gridpoint that will be
     changed to land.
  2. Define cloning gridpoint. Attributes from this point will be copied to the water body points.
     This point should be representative of the land cover that would be found where the body is
     located.
  3. Execute code. The output filename should be changed to reflect the WPS geogrid configuration.
  4. Place output file in WRF/WPS directory.
  5. Run metgrid.exe, real.exe, and wrf.exe. The model should run without the water body present. 

Acknowledgements

Resources for this project were donated by the Department of Geography and Meteorology at Valparaiso University

Important information about what fields needed to be modified came from the following website: http://www.wrfems.info/viewtopic.php?t=141 

Code for copying netCDF variables and dimensions from one file to another in Python came from the following website: http://stackoverflow.com/questions/15141563/python-netcdf-making-a-copy-of-all-variables-and-attributes-but-one 
