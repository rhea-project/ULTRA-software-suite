## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import math as m
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Anti_icing system Weight Calculation
# ----------------------------------------------------------------------

def weight_anti_icing (SPAN ,SWEEP,FNAC,FNENG,WF):
    """ Calculate the weight of anti_icing system based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
       SPAN    - Wing span                                               [m]
       SWEEP  - Quarter chord sweep angle of wing                       [deg]
       FNAC   - Average diameter of each scaled engine                  [m]
       FNENG  - Total number of engines
       WF     - Maximum fuselage width                                  [m]
    Outputs:
        WAI      - Weight of anti_icing system                          [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    SPAN =SPAN /Units.ft 
    FNAC=FNAC/Units.ft 
    WF=WF/Units.ft   
    SWEEP=m.radians(SWEEP)
    print(SWEEP)
    WAI=(SPAN/m.cos(SWEEP))+3.8*FNAC*FNENG+1.5*WF
    WAI=WAI*Units.lb  
    return WAI
