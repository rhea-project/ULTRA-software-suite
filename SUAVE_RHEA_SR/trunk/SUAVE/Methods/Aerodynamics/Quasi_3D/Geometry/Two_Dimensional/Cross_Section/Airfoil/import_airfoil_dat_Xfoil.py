## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil
# import_airfoil_dat_Xfoil.py
# 
# Created:  
# Modified: Sep 2016, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np

# ------------------------------------------------------------
#  import airfoil dat
# ------------------------------------------------------------

## @ingroup Methods-Geometry-Two_Dimensional-Cross_Section-Airfoil
def import_airfoil_dat_Xfoil(filename):
    """Import an airfoil data file and stores it in a numpy array.
    
    Assumptions:
    Airfoil file in Xfoil format

    Source:
    None

    Inputs:
    filename   <string>

    Outputs:
    data       numpy array with airfoil data

    Properties Used:
    N/A
    """     
    
    filein = open(filename,'r')
    data = {}   
    data['header'] = filein.readline().strip() 
    
    sections = ['airfoil']
    data['airfoil'] = []
    
    while True:
        line = filein.readline()
        if not line: 
            break
        
        line = line.strip()
        
        point = list(map(float,line.split()))
        data['airfoil'].append(point)
        
    for k,v in data.items():
        data[k] = np.array(v)
        
    return data


