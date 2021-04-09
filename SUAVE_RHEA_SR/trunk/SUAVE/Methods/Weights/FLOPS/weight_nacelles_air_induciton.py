## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import math as m
import pandas as pd
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Nacelles/Air induction Weight Calculation
# ----------------------------------------------------------------------

def weight_nac(ac_type,TNAC,thrust,thrso,NEF,WF,DF,XNAC,DNAC):
    
    """ Calculate nacelles/air induction weight of the aircraft based on the FLOPS methods
          
    Inputs:
        ac_type :  aircraft type: 1 for transport aircraft
                                  2 for fighter/attack aircrat
        TNAC    : number of nacelles                          
        thrust  : rated thrust for each scaled engine            [N]   
        thrso   : rated thrust for each baseline engine          [N]
        NEF     : number of fuselage-mounted engines             [Unitless]
        WF      : max fuselage width                             [m]
        DF      : max fuselage depth                             [m]
        
    Outputs:
       weight_nac                                                [kg]      
        
    Properties Used:
        N/A
    """ 
    
    # unpack inputs
    XNAC   = XNAC * 3.28084                  #m to ft
    DNAC   = DNAC * 3.28084                 #m to ft   
    thrust = thrust*0.22481                 #N to lb
    thrso  = thrso*0.22481                  #N to lb

    if XNAC == 0:
        XNAC       = 0.07*m.sqrt(thrso)
    if DNAC == 0:
        DNAC       = 0.04*m.sqrt(thrso)

    # for transport aircraft
    if ac_type=='transport' or ac_type=='general aviation':
        weight_nac = 0.25*TNAC*XNAC*DNAC*thrust**0.36
    else:
        weight_nac = 0.25*TNAC*XNAC*DNAC*thrust**0.36 + \
                     1.06*(thrust*NEF)**0.23*(WF+DF)**1.4*Mmax**0.83
        
    weight_nac = weight_nac*Units.lb
    
    return weight_nac
