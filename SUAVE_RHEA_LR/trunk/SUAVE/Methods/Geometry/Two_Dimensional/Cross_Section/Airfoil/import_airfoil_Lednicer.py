## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil
# import_airfoil_dat.py
# 
# Created:  Jun 2020 S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np

# ------------------------------------------------------------
#  import airfoil dat
# ------------------------------------------------------------

## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil
def import_airfoil_Lednicer(file_name):
    """Import an airfoil data file in the Lednicer format and stores it in a numpy array.
    
    Assumptions:
    Airfoil file in Lednicer format

    Source:
    None

    Inputs:
    filename   <string>

    Outputs:
    data       numpy array with airfoil data

    Properties Used:
    N/A
    """     
    
    airfoil_coord_file = open(file_name, "r")

    # Read the header
    header = airfoil_coord_file.readline().replace("\n"," ").replace("\t"," ").split()

    # Read the body
    body = airfoil_coord_file.read().replace("\n"," ").replace("\t"," ").split()
    w = 0
    points = int(len(body)/2)
    x1 = np.zeros(points)
    y1 = np.zeros(points)
    airfoil_coord_file.close()
    airfoil_coord_file = open(file_name, "r")
    header = airfoil_coord_file.readline().replace("\n"," ").replace("\t"," ").split()

    for i in range(points):
        coordinate_point = airfoil_coord_file.readline().replace("\n"," ").replace("\t"," ").split()
        x1[i] = coordinate_point[0]
        y1[i] = coordinate_point[1]

    return x1, y1, points
