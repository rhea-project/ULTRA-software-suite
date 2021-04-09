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
#   Air Conditioning Weight Calculation
# ----------------------------------------------------------------------

def weight_air_cond(ac_type,FPAREA,DF,NPASS,VMAX,WAVONC,FENG,FTHRST):
    """ Calculate the weight of air conditioning system based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                     # 1 for transport/general aviation aircraft
                                                                     # 2 for fighter/attack aircraft
        FPAREA   - Fuselage planform area                           [m²]
        DF       - Maximum fuselage depth                           [m]
        NPASS    - Total number of passengers                                                            
        VMAX     - Maximum Mach number
        WAVONC   - Weight of the avionics system group              [kg]
        FENG     - Total number of engines
        FTHRST   - Rated thrust of each scaled engine               [kN]
    Outputs:
        WAC      - Weight of the air conditioning system            [kg]
    Properties Used: 
        N/A
    """ 
     # unpack inputs
    FPAREA=FPAREA/Units.ft**2
    DF =DF /Units.ft
    WAVONC=WAVONC/Units.lb
    FTHRST=FTHRST*224.81  #kN to lb

    if ac_type==1:
        WAC=(3.2*(FPAREA*DF)**0.6+9*NPASS**0.83)*VMAX+0.075*WAVONC
    elif ac_type==2:
        WAC=0.075*WAVONC+0.37*FENG*FTHRST**0.6*VMAX**0.57
    WAC=WAC*Units.lb
    return WAC
