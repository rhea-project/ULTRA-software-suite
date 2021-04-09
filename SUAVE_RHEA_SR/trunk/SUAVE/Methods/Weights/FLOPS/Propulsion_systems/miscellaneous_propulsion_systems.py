## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S.Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import math as m
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#  Engine Weight Calculation
# ----------------------------------------------------------------------

def weight_engine_miscellaneous(NEW,NEF,ac_type,THRUST,FNAC,WPMISC,NFCLR,VMAX):

    """ Calculate the miscellaneous propulsion systems weight  based on the FLOPS methods    
   
    Source: 
        N/A
        
    Inputs:
    
        NEW     - number of wing mounted engines
        NEF     - number of fuselage mounted engines.                                    The default value is 0 .
        ac_type - 1 for transport aircraft, 2 for fighter
        THRUST  - Rated thrust of each scaled engine                               [N]
        FNAC    - Average diameter of each scaled engine                           [m]
        WPMISC  - Additional miscellaneous propulsion system weight               [kg]
        NFCLR   - number of flight crew. For fighter equal 1.
        VMAX    - Maximum Mach number
        
    Outputs:
        WPMSC = Total miscellaneous propulsion systems weight                     [kg]                                             
      
    Properties Used:
        N/A
    """ 
    # unpack inputs
    THRUST = THRUST*0.22481               #N to lb
    FNAC   = FNAC/Units.ft                #m to ft
    WPMISC = WPMISC/Units.lb
    NENG   = NEW+NEF
    
    #Calculate the engine controls weight
    if ac_type=='transport' or ac_type=='general aviation':
       WEC    = 0.26*NENG*THRUST**0.5
       WSTART = 11*NENG*VMAX**0.32*FNAC**1.6
    elif ac_type=='fighter':
       WEC    = 0.106*(NENG*THRUST*NFCLR)**0.55  
       WSTART = 0.0072*THRUST*NENG
       
    WPMSC = WEC+WSTART+WPMISC
    WPMSC = WPMSC*Units.lb
    
    return WPMSC

