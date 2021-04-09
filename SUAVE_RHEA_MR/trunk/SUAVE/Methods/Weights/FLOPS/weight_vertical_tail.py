## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import math as m
from SUAVE.Core import Units
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Vertical Tail Weight Calculation
# ----------------------------------------------------------------------

def weight_vt(DG,TRVT,NVERT,SVT,ULF,ARVT,SWPVT,HHT,TCVT,p_h,p_0,VCMN,AC_type):

    """ Calculate the vertical tail weight of the aircraft based on the FLOPS methods
            
    Inputs:
        DG      - Design gross weight                             [kg] 
        TRVT    - Vertical tail theoretical taper ratio
        NVERT   - Number of vertical tails
        SVT     - Vertical tail theoretical area                  [m**2]
        ULF     - Structural ultimate load factor 
        ARVT    - Vertical tail theoretical aspect ratio
        SWPVT   - Vertical tail sweep angle at 25% chord
        HHT     - Horizontal tail mounting location indicator #ranging from 0.0 for body mounted to 1.0 for a T-tail
        TCVT    - Vertical tail thickness to chord ratio  #assumed equal to 10% for all aircrafts
        p_h     - pressure at cruise altitude                    [Pa]
        p_0     - pressure at sea level                          [Pa]
        VCMN    - Cruise Mach number
        AC_type - Aircraft type   

    Outputs:
       WVT                                                       [kg]
        
    Properties Used:
        N/A
    """ 
    
    # unpack inputs
    DG  = DG/Units.lb           # Convert kg to lb 
    SVT = SVT/Units.ft**2       # Convert meters squared to ft squared
   
    #Calculate weight of the horizontal tail using  FLOPS methods
    #Transport aircraft
    if AC_type == 'transport': 
      transport_ac_WVT = 0.32*DG**0.3*(TRVT+0.5)*NVERT**0.7*SVT**0.85 
      WVT              = transport_ac_WVT
     #fighter/attack aircraft 
    elif AC_type == 'fighter':
        fighter_ac_WVT = 0.212*DG**0.3*(TRVT+0.5)*NVERT**0.7*SVT**0.94*ARVT**0.5/((m.cos(SWPVT))**1.5)
        WVT            = fighter_ac_WVT
    #general aviation aircraft 
    else:
         DELTA                = p_h/p_0
         QCRUS                = 1481.35*DELTA*VCMN**2
         CSVT                 = m.cos(SWPVT)
         general_aviation_WVT = 0.073*(1+0.2*HHT)*(ULF*DG)**0.376*QCRUS**0.122*SVT**0.873*((ARVT/(CSVT**2))**0.357)/((100*TCVT)/CSVT)**0.49
         WVT                  = general_aviation_WVT
         
    WVT = WVT*Units.lb          # Convert lb to kg
    
    return WVT





