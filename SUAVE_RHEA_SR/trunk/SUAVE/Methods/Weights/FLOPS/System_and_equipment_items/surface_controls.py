## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S.Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


import pandas as pd
from SUAVE.Core import Units

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Flight controls Weight Calculation
# ----------------------------------------------------------------------

def weight_surface_controls(ac_type,DG,VMAX,SFLAP,p_0,p_h,ULF,SW):
    
    """ Calculate the surface controls weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type - aircraft type                
        VMAX    - maximum Mach number
        SFLAP   - total movable wing surface area including flaps, elevators, spoilers, etc. [m²]
        DG      - Design gross weight                                                        [kg] 
        p_0     - Pressure at sea level                                                      [Pa]
        p_h     - Pressure at flight altitude H                                              [Pa]
        ULF     - Structural ultimate load                                                            #default value is equal to 3.75
        SW      - Wing reference area                                                        [m²]

    Outputs:
      WSC       - weight of the surface controls               [kg]
      
    Properties Used:
        N/A
    """ 
    # unpack inputs

    DG    = DG/Units.lb                     # Convert kg to lb 
    SFLAP = SFLAP/Units.ft**2               # Convert meters squared to ft squared    
    SW    = SW/Units.ft**2                  # Convert meters squared to ft squared

    #Calculate the surface controls weight using FLOPS methods
    if ac_type=='transport':
        WSC = 1.1*(VMAX**0.52)*SFLAP**0.6*DG**0.32
    elif ac_type=='fighter':
        WSC = 2.95*SFLAP**0.45*DG**0.36
    else:
        delta = p_h/p_0
        QDIVE = 1481.35*delta*VMAX**2
        WSC   = 0.404*SW**0.317*(DG/1000)**0.602*ULF**0.525*QDIVE**0.345

    WSC = WSC*Units.lb
    
    return WSC

