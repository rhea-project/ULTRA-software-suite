## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
from SUAVE.Core import Units

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Canard Weight Calculation
# ----------------------------------------------------------------------

def weight_canard(DG,can_area,can_TR):
    """ Calculate the canard weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        DG       - Design gross weight                      [kg] 
        can_area - Fin theoretical area                     [kg]
        can_TR   - Fin theoretical taper ratio
       
    Outputs:
      weight_fin                                                     
        
    Properties Used:
        N/A
    """ 
     # unpack inputs
    can_area = can_area/Units.ft**2 # Convert meters squared to ft squared
    DG       = DG/Units.lb          # Convert kg to lb

    #Calculate the canard weight using FLOPS methods
    weight_can = 0.53*can_area*DG**0.2*(can_TR+0.5)
    weight_can = weight_can*Units.lb
    
    return weight_can
