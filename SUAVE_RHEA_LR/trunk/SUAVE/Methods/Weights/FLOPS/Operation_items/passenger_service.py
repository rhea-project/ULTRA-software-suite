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
#   Passenger service Weight Calculation
# ----------------------------------------------------------------------

def weight_pass_service(NPF,NPB,NPT,design_range,VMAX):
    
    """ Calculate the weight reserved for the passenger service  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
      NPF          - Number of first class passengers
      NPB          - Number of business class passengers
      NPT          - Number of tourist class passengers
      design_range - Desing range                                                [m]
      VMAX         - Maximum Mach number
       
ns
    Outputs:
      WSRV - Weight of passenger service for transport aircraft                 [kg]
      
    Properties Used:
        N/A
    """ 
    
    # unpack inputs
    design_range = 0.000539957*design_range
     
    WSRV = (5.164*NPF+3.846*NPB+2.529*NPT)*((design_range/VMAX)**0.225)
    WSRV = WSRV*Units.lb
    
    return WSRV
