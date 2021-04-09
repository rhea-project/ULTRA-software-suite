## @ingroup Methods-Aerodynamics-Xfoil
#write_geometry.py
# 
# Created:  July 2019, S. Karpuk
# Modified:
# 


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
from SUAVE.Core import Units, Data
from .purge_files import purge_files
from SUAVE.Methods.Aerodynamics.XFOIL.Data.Settings    import Settings
import numpy as np


## @ingroup Methods-Aerodynamics-AVL
def write_geometry(xfoil_object):
    """This function writes the translated airfoil geometry into text file read 
    by Xfoil when it is called

    Assumptions:
        None
        
    Source:
        Drela M., Xfoil, https://web.mit.edu/drela/Public/web/xfoil

    Inputs:
        xfoil_object

    Outputs:
        None

    Properties Used:
        N/A
    """

    
    # Configure the Data structure settings
    xfoil_object.settings = Settings()

    # unpack inputs
    n_segments = len(xfoil_object.airfoils.Segments.keys())
    for i in range(n_segments):
        with open(xfoil_object.airfoils.Segments[i].Airfoil.airfoil.coordiante_file,'w') as airfoil_file:
            airfoil_file.write(xfoil_object.airfoils.Segments[i].tag + '\n')
            np.savetxt(airfoil_file, xfoil_object.airfoils.Segments[i].Airfoil.airfoil.points, fmt="%6.4f")

    return
