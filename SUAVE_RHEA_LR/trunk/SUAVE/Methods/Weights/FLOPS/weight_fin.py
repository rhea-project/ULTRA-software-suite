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
#   Fin Weight Calculation
# ----------------------------------------------------------------------

def weight_fin(DG,fin_area,fin_TR,n_fin):
    """ Calculate the fin weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        DG       - Design gross weight                      [kg] 
        fin_area - Fin theoretical area                     [kg]
        fin_TR   - Fin theoretical taper ratio
        n_fin    - Number of fins

    Outputs:
      weight_fin                                            [kg]        
        
    Properties Used:
        N/A
    """ 
     # unpack inputs
    fin_area = fin_area/Units.ft**2 # Convert meters squared to ft squared
    DG       = DG/Units.lb          # Convert kg to lb


    #Calculate the fin weight using FLOPS methods
    weight_fin = 0.32*DG**0.3*fin_area**0.85*(fin_TR+0.5)*n_fin
    weight_fin = weight_fin*Units.lb
    return weight_fin 
