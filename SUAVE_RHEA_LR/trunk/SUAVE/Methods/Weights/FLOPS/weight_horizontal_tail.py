## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Horizontal Tail Weight Calculation
# ----------------------------------------------------------------------

def weight_ht(SHT,DG,TRHT,ULF,VCMN,p_h,p_0,AC_type):
    
    """ Calculate the horizontal tail weight of the aircraft based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        SHT  -Horizontal tail theoretical area                [m**2]
        DG   -Design gross weight                             [kg] 
        TRHT -Horizontal tail theoretical taper ratio
        ULF  -Structural ultimate load factor                 #  Default value 2.5
        p_h   -pressure at cruise altitude                   [Pa]
        p_0   -pressure at sea level                         [Pa]
        VCMN -Cruise Mach number         
        AC_type -Aircraft type  

    Outputs:
       WHT                                                    [kg]
        
    Properties Used:
        N/A
    """ 
     # unpack inputs
    SHT = SHT/Units.ft**2                       # Convert meters squared to ft squared
    DG  = DG/Units.lb                           # Convert kg to lb

    #Calculate weight of the horizontal tail FLOPS methods
    #Transport aircraft
    if AC_type=='transport': 
      transport_ac_WHT = 0.53*SHT*DG**0.2*(TRHT+0.5) 
      WHT              = transport_ac_WHT
    #fighter/attack aircraft 
    elif AC_type=='fighter':
        fighter_ac_WHT = 0.002*SHT**0.87*(ULF*DG)**0.66
        WHT            = fighter_ac_WHT
    #general aviation aircraft 
    else:
         DELTA                = p_h/p_0
         QCRUS                = 1481.35*DELTA*VCMN**2
         general_aviation_WHT = 0.016*SHT**0.873*(ULF*DG)**0.414*QCRUS**0.122
         WHT                  = general_aviation_WHT
         
    WHT = WHT*Units.lb                          # Convert lb to kg


    return WHT
