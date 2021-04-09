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
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Instruments Weight Calculation
# ----------------------------------------------------------------------

def weight_instruments(ac_type,FPAREA,VMAX,NFLCR,FNEW,FNEF,XF,DF,NFUSE):
    
    """ Calculate the instrument weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Airfract type                                    
        FPAREA   - Fuselage planform area                            [m²]
        VMAX     - Maximum Mach number
        NFLCR    - Number of flight crew
        FNEW     - Number of wing mounted engines
        FNEF     - Number of fuselage mounted engines     
        XF       - Total fuselage length                             [m]
        DF       - Maximum fuselage depth                            [m]                     
        NFUSE    - Number of fuselages

    Outputs:
      WIN      - Weight of the instrumenmt system group             [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    FPAREA = FPAREA/Units.ft**2 
    XF     = XF/Units.ft                  # Convert meter to ft
    DF     = DF/Units.ft                  # Convert meter to ft

    FPAREA = FPAREA * NFUSE

    #Calculate the instrument weight using FLOPS methods
    if ac_type=='transport':
        WIN = 0.48*FPAREA**0.57*VMAX**0.5*(10+2.5*NFLCR+FNEW+1.5*FNEF)
    elif ac_type=='fighter':
      WIN = 0.09*NFUSE*XF*DF*(1+2.5*NFLCR+0.1*FNEW+0.15*FNEF)
      
    WIN = WIN*Units.lb
    
    return WIN
