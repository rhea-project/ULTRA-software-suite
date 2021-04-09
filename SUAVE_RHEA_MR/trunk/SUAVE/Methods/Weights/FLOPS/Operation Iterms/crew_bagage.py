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
#   Crew and baggage Weight Calculation
# ----------------------------------------------------------------------

def weight_crew_baggage(NPASS,
                        ac_type,
                        CARBAS=0):
    """ Calculate the crew and baggage weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        NPASS - number of passenger
        ac_type  - aircraft type                 : 1 for Transport/HWB
                                                   2 for fighter/attack aircraft
        CARBAS   - Carrier based aircraft switch : 1 for carrier-based aircraft
                                                   0 for land-based aircraft
       
ns
    Outputs:
      WSTUAB+WFLCRB    - Total  Weight of flight attendants ,galley crew , baggage and flight crew  [kg]                                             
      
    Properties Used:
        N/A
    """ 
    

    #Calculate the crew and baggage weight using FLOPS methods
    if NPASS<0 :
     NSTU=0
    elif NPASS<51:
     NSTU=1
    else:
     NSTU=1+ m.ceil(NPASS/40)
    NGALC=0 if NPASS<151 else (1+ m.ceil(NPASS/250))
    if (ac_type==1):
      NFLCR=2 if  NPASS<151 else 3
    else:
     NFLCR=1
    WSTUAB=NSTU*155+NGALC*200
    WFLCRB= NFLCR*(215-35*CARBAS) if ac_type==2 else NFLCR*(225-35*CARBAS)
    WSTUAB=WSTUAB*Units.lb
    WFLCRB=WFLCRB*Units.lb
    return WSTUAB+WFLCRB
    
      
