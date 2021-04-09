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
#   Avionics Weight Calculation
# ----------------------------------------------------------------------

def weight_avionics(ac_type,DESRNG,NFLCR,FPAREA,NFUSE,XL,DF,VMAX,CARBAS):
    
    """ Calculate the weight of avionics system based on the FLOPS methods
       
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                   
        DESRNG   - Design range                                     [m]
        NFLCR    - Number of flight crew
        FPAREA   - Fuselage planform area                           [m²]
        NFUSE    - Number of fuselages
        XL       - Total fuselage length                            [m]
        DF       - Maximum fuselage depth                           [m]
        CARBAS   - Carrier based aircraft switch                     # 1 for carrier-based aircraft
                                                                     # 0 for land-based aircraft
        VMAX     - Maximum Mach number
        
    Outputs:
      WAVONC     - Weight of the instrumenmt system group             [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    DESRNG = 0.000539957 * DESRNG
    FPAREA = FPAREA/Units.ft**2
    XL     = XL/Units.ft
    DF     = DF/Units.ft 

    FPAREA = FPAREA * NFUSE

    if ac_type=='transport':
        WAVONC = 15.8*DESRNG**0.1*NFLCR**0.7*FPAREA**0.43
    elif ac_type=='fighter':
        WAVONC = 0.41*(NFUSE*XL*DF)**1.3*(1+0.37*CARBAS)*VMAX
        
    WAVONC = WAVONC*Units.lb

    return WAVONC
    
