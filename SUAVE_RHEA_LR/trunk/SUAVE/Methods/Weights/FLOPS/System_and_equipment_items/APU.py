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
#   Auxiliary Power Unit Weight Calculation
# ----------------------------------------------------------------------

def weight_APU(NPASS,FPAREA,NFUSE):
    
    """ Calculate the  weight of the APU  based on the FLOPS methods
       
    Source: 
        N/A
        
    Inputs:       
        
        FPAREA  - Fueselage planform area                           [m²]
        NPASS   - Total number of passengers
        NFUSE   - Number of fuselages

    Outputs:
      WAPU      - Weight of the APU                                 [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    FPAREA  = FPAREA /Units.ft**2  # Convert meter to inches

    FPAREA  = FPAREA * NFUSE
    
    #Calculate the APU weight using FLOPS methods
    WAPU = 54*FPAREA**0.3+5.4*NPASS**0.9
    WAPU = WAPU*Units.lb
    
    return WAPU


