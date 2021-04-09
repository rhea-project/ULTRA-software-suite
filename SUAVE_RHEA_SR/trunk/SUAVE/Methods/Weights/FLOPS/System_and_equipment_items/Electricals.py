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
#   Electrical system Weight Calculation
# ----------------------------------------------------------------------

def weight_electrical(ac_type,XL,WF,NFUSE,FNENG,NFLCR,NPASS,VMAX,SPAN):
    
    """ Calculate the weight of the electrical system based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                 
        XL       - total fueslage length                            [m]
        WF       - Maximum fuselage width                           [m]
        NFUSE    - Number of fuselages
        FNENG    - Total number of engines
        NFLCR    - Number of flight crew
        NPASS    - Total number of passengers
        VMAX     - Maximum Mach number
        SPAN     - Wing span                                           #used for the fighter/attack aircraft
    Outputs:
      WIN      - Weight of the instrumenmt system group             [kg]
    Properties Used: 
        N/A
    """
    
    # unpack inputs
    XL   = XL/Units.ft               # Convert meter to ft
    WF   = WF/Units.ft               # Convert kg to lb
    SPAN = SPAN/Units.ft

    if ac_type=='transport' or ac_type=='general_aviation':
        WELEC = 92*XL**0.4*WF**0.14*NFUSE**0.27*FNENG**0.69*(1+0.044*NFLCR+0.0015*NPASS)
    elif ac_type=='fighter':
        WELEC = 10*(XL+SPAN)**0.85*NFUSE**0.27*VMAX**0.1*(1+0.1*NFLCR)
        
    WELEC = WELEC*Units.lb
    
    return WELEC
