## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

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

def weight_avionics(ac_type,
                    DESRNG,
                    NFLCR,
                    FPAREA,
                    NFUSE,
                    XL,
                    DF,
                    CARBAS,
                    VMAX):
    """ Calculate the weight of avionics system based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                     # 1 for transport/general aviation aircraft
                                                                     # 2 for fighter/attack aircraft
        DESRNG   - Design range                                     [nmi]
        NFLCR    - Number of flight crew
        FPAREA   - Fuselage planform area                           [m²]
        NFUSE    - Number of fuselages
        XL       - Total fuselage length                            [m]
        DF       - Maximum fuselage depth                           [m]
        CARBAS   -Carrier based aircraft switch                     # 1 for carrier-based aircraft
                                                                    # 0 for land-based aircraft
        VMAX     - Maximum Mach number
        
    Outputs:
      WAVONC     - Weight of the instrumenmt system group             [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    FPAREA=FPAREA/Units.ft**2
    XL=XL/Units.ft
    DF=DF/Units.ft 

    if ac_type==1:
        WAVONC=15.8*DESRNG**0.1*NFLCR**0.7*FPAREA**0.43
    elif ac_type==2:
        WAVONC=0.41*(NFUSE*XL*DF)**1.3*(1+0.37*CARBAS)*VMAX
    WAVONC=WAVONC*Units.lb  
    return WAVONC
    