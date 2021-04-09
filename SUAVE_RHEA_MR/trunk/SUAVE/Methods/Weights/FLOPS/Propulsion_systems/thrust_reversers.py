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

def weight_thrust_reversers(NEW,NEF,THRUST):
    
    """ Calculate the thrust reversers weight  based on the FLOPS methods
       
    Source: 
        N/A
        
    Inputs:
    
        NEW    - number of wing mounted engines
        NEF    - number of fuselage mounted engines.                                The default value is 0 .
        THRUST - Rated thrust of each scaled engine                       [N]
        
    Outputs:
        WTHR = Weight of the thrust reversers                             [kg]                                             
      
    Properties Used:
        N/A
    """ 
    # unpack inputs
    THRUST = THRUST*0.22481                     #N to lb
    NENG   = NEW+NEF
    
    #Calculate the total number of nacelles
    TNAC = NENG+0.5*(NENG-2*m.floor(NENG/2))
    WTHR = 0.034*THRUST*TNAC
    WTHR = WTHR*Units.lb
    
    return WTHR

