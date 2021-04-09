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
#  engine oill Weight Calculation
# ----------------------------------------------------------------------

def weight_engine_oil(rat_thrust,num_eng):
    
    """ Calculate the engine oil weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
    
        num_eng - number of engines                                         [Unitless]
        
        rat_thrust - Rated thrust of each scaled engine                     [N]
        
    Outputs:
      EOW - Weight of the engine oil                                        [kg]                                             
      
    Properties Used:
        N/A
    """
    
    # unpack inputs
    rat_thrust=rat_thrust*0.22481 #kN to lb

    #Calculate the unusable fuel  weight using FLOPS methods
    EOW = 0.082*num_eng*rat_thrust**0.65
    EOW = EOW*Units.lb
    
    return EOW


