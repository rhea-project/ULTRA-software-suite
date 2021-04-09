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
#   Air Conditioning Weight Calculation
# ----------------------------------------------------------------------

def weight_air_cond(ac_type,FPAREA,DF,NPASS,VMAX,WAVONC,FENG,FTHRST,NFUSE):
    
    """ Calculate the weight of air conditioning system based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                     
        FPAREA   - Fuselage planform area                           [m²]
        DF       - Maximum fuselage depth                           [m]
        NPASS    - Total number of passengers                                                            
        VMAX     - Maximum Mach number
        WAVONC   - Weight of the avionics system group              [kg]
        FENG     - Total number of engines
        FTHRST   - Rated thrust of each scaled engine               [kN]
        NFUSE    - Number of fuselages
    Outputs:
        WAC      - Weight of the air conditioning system            [kg]
    Properties Used: 
        N/A
    """ 
     # unpack inputs
    FPAREA = FPAREA/Units.ft**2
    DF     = DF/Units.ft
    WAVONC = WAVONC/Units.lb
    FTHRST = FTHRST*0.22481          #N to lb

    FPAREA = FPAREA * NFUSE

    if ac_type=='transport':
        WAC = (3.2*(FPAREA*DF)**0.6+9*NPASS**0.83)*VMAX+0.075*WAVONC
    elif ac_type=='fighter':
        WAC = 0.075*WAVONC+0.37*FENG*FTHRST**0.6*VMAX**0.57
        
    WAC = WAC*Units.lb
    
    return WAC
