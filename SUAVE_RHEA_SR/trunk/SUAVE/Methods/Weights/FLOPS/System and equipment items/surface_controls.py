## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


import pandas as pd

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Flight controls Weight Calculation
# ----------------------------------------------------------------------

def weight_surface_controls(ac_type,DG,VMAX,SFLAP,p_0,p_h,ULF,SW):
    """ Calculate the surface controls weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type - aircraft type                 : 1 for Transport aircraft
                                                   2 for fighter/attack aircraft
                                                   3 for general aviation aircraft
        VMAX    - maximum Mach number
        SFLAP   - total movable wing surface area including flaps, elevators, spoilers, etc. [ft²]
        DG      - Design gross weight                                       [lb] 
        p_0     - Pressure at sea level                                     [psi]
        p_h     - Pressure at flight altitude H                             [psi]
        ULF     - Structural ultimate load        #default value is equal to 3.75
        SW      - Wing reference area                                       [ft²]
ns
    Outputs:
      WSC       - weight of the surface controls               [lb]
    Properties Used:
        N/A
    """ 
    

    #Calculate the surface controls weight using FLOPS methods
    if ac_type==1:
        WSC=1.1*(VMAX**0.52)*SFLAP**0.6*DG**0.32
    elif ac_type==2:
        WSC=2.95*SFLAP**0.45*DG**0.36
    else:
     delta=p_h/p_0
     QDIVE= 1481.35*delta*VMAX**2
     WSC=0.404*SW**0.317*(DG/1000)**0.602*ULF**0.525*QDIVE**0.345
    return WSC

