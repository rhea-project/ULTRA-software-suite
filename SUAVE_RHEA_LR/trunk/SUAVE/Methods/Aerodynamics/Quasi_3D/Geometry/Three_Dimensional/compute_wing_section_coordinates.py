## @ingroup Methods-Geometry-Three_Dimensional
# compute_wing_section_coordinates.py
# 
# Created:  Jul 2019, S. Karpuk, 
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data
from SUAVE.Methods.Geometry.Two_Dimensional.Cross_Section.Airfoil import import_airfoil_dat_Xfoil
from SUAVE.Methods.Aerodynamics.Quasi_3D.transform_section import transform_section
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
#  Compute Wing Section Coordinates
# ----------------------------------------------------------------------

## @ingroup Methods-Geometry-Three_Dimensional
def compute_wing_section_coordinates(avl_object,wing):
    """Create a numpy array of wing sections presented in the aircraft.
    
    Assumptions:

    Source:
    None

    Inputs:
    wing

    Outputs:

    Properties Used:
    N/A
    """

    # create designed airfoils with similar number of points    
    if len(wing.Segments.keys())>0:

        n_segments  = len(wing.Segments.keys())
        
        for i_segs in range(n_segments):
            if i_segs == 0:
                root_airfoil = import_airfoil_dat_Xfoil(wing.Segments[i_segs].Airfoil.airfoil.coordinate_file)
                wing.Segments[i_segs].Airfoil.airfoil.points = root_airfoil['airfoil']
            else:
                section_airfoil  = import_airfoil_dat_Xfoil(wing.Segments[i_segs].Airfoil.airfoil.coordinate_file)
                airfoil_points = interpolate_airfoil_coordinates(root_airfoil, section_airfoil)
                wing.Segments[i_segs].Airfoil.airfoil.points = airfoil_points

        transform_section(avl_object,wing)

    return


def interpolate_airfoil_coordinates(root_airfoil, section_airfoil):
    """Create a numpy array of wing sections with similar point distribution as in the root airfoil.
    
    Assumptions:
    Airfoil file in Xfoil format

    Source:
    None

    Inputs:
    root_airfoil             numpy array with airfoil data
    section_airfoil          numpy array with airfoil data

    Outputs:
    data                     numpy array with airfoil data

    Properties Used:
    N/A
    """    

    fig, ax = plt.subplots()
    x1 = root_airfoil['airfoil'][:,0]
    y1 =root_airfoil['airfoil'][:,1]
    x2 = section_airfoil['airfoil'][:,0]
    y2 = section_airfoil['airfoil'][:,1]

    x2u = []
    y2u = []
    x2l = []
    y2l = []
    x1u = []
    x1l = []

    # split the surface of the root airfoil
    x1u,x1l,y1u,y1l = split_airfoil(x1,y1)
    
    # split the surface of the section airfoil
    x2u,x2l,y2u,y2l = split_airfoil(x2,y2)
 
    # interpolate the airfoils
    y3u  = np.interp(x1u[::-1],x2u[::-1],y2u[::-1])
    y3l  = np.interp(x1l,x2l,y2l)
    y3   = np.concatenate((y3u[::-1], y3l), axis=None)
    x11  = np.array([x1])
    y31  = np.array([y3])
    data = np.concatenate((x11.T, y31.T), axis=1)

    return data

def split_airfoil(x,y):

    xu = []
    xl = []
    yu = []
    yl = []
    
    i = 0
    while x[i] != 0:
        xu.append(x[i])
        yu.append(y[i])
        i = i + 1
    for j in range(i,np.shape(x)[0]):
        xl.append(x[j])
        yl.append(y[j])

    return  xu,xl,yu,yl

