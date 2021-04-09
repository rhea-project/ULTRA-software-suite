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
#   Furnishings and Equipment Weight Calculation
# ----------------------------------------------------------------------

def weight_furnishings(ac_type,
                       NFLCR,
                       NPF,
                       NPB,
                       NPT,
                       XLP,
                       NFUSE,
                       VMAX,
                       XL,
                       CARGF,
                       ACABIN,
                       WF,
                       DF,
                       NBAY,
                       SWPLE):
    """ Calculate the weight of the furnishings and equipment based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        ac_type  - Aircraft type                                     # 1 for transport/general aviation aircraft
                                                                     # 2 for fighter/attack aircraft
                                                                     # 3 for HWB aircraft
        NFLCR    - Number of flight crew
        NPF      - Number of first class passengers
        NPB      - Number of business class passengers
        NPT      - Number of tourist class passengers
        XLP      - Length of passenger compartmenet                 [m]
        NFUSE    - Number of fuselages
        VMAX     - Maximum Mach number
        XL       - Total fuselage length                            [m] 
        CARGF    - Cargo aircraft floor factor                     # used for HWB aircraft 0 for passenger transport 
                                                                   #                      1 for military cargo transport
        ACABIN   - Passenger cabin floor area                       [m]
        WF       - Maximum fuselage width                           [m]
        DF       - Maximum fuselage depth                           [m]
        NBAY     - Number of passenger bays
        SWPLE    - Sweep angle of the passenger cabin               [deg]
    Outputs:
      WFURN      - Weight of the furnishings group                  [kg]
    Properties Used: 
        N/A
    """ 
     # unpack inputs
    XLP=XLP/Units.ft  # Convert meter to ft
    XL=XL/Units.ft   # Convert kg to lb
    ACABIN=ACABIN/Units.ft**2
    WF=WF/Units.ft 
    DF=DF/Units.ft 
    if ac_type==1:
        WFURN=127*NFLCR+112*NPF+78*NPB+44*NPT+2.6*XLP*(WF+DF)*NFUSE
    elif ac_type==2:
        WFURN=80*NFLCR*VMAX**0.38*XL**0.25
    else:
        WFURN=127*NFLCR+112*NPF+78*NPB+44*NPT+2.6*(1-CARGF)*((ACABIN*(WF+DF*NBAY)/WF)+WF*DF*(1+(1/m.cos(SWPLE))))
    WFURN=WFURN*Units.lb    
    return WFURN

