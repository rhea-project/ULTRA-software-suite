## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


import pandas as pd
import math as m
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Crew and baggage Weight Calculation
# ----------------------------------------------------------------------

def weight_crew_baggage(NPASS,ac_type,CARBAS=0):
    
    """ Calculate the crew and baggage weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        NPASS - number of passenger                                              [Unitless]
        ac_type  - aircraft type                                                 [Unitless]                 

        CARBAS   - Carrier based aircraft switch : 1 for carrier-based aircraft
                                                   0 for land-based aircraft
       
ns
    Outputs:
      WSTUAB+WFLCRB    - Total  Weight of flight attendants ,galley crew , baggage and flight crew  [kg]                                             
      
    Properties Used:
        N/A
    """ 
    

    #Calculate the crew and baggage weight using FLOPS methods
    if NPASS < 0 :
     NSTU = 0
    elif NPASS < 51:
     NSTU = 1
    else:
     NSTU = 1+ m.ceil(NPASS/40)
     
    if NPASS<151: 
        NGALC = 0
    else:
        NGALC = (1+ m.ceil(NPASS/250))
        
    if ac_type=='transport' or ac_type=='HWB':
        if  NPASS < 151:
            NFLCR = 2
        else:
            NFLCR = 3
    else:
        NFLCR = 1
        
    WSTUAB = NSTU*155+NGALC*200
    if ac_type=='fighter':
        WFLCRB = NFLCR*(215-35*CARBAS)
    else:
        WFLCRB = NFLCR*(225-35*CARBAS)
        
    WSTUAB = WSTUAB*Units.lb
    WFLCRB = WFLCRB*Units.lb
    
    return WSTUAB+WFLCRB
    
      
