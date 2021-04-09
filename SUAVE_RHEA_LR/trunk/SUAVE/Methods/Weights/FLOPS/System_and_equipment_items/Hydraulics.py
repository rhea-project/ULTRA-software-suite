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
#   Hydraulic system Weight Calculation
# ----------------------------------------------------------------------

def weight_hydraulics(FPAREA,VMAX,SW,FNEW,FNEF,VARSWP,NFUSE,HYDPR=3000):
    
    """ Calculate the weight of the hydraulic system  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
       
        FPAREA   - Fuselage planform area                           [m²]
        VMAX     - Maximum Mach number
        SW       - Wing reference area                              [m²]
        FNEW     - Number of wing mounted engines
        FNEF     - Number of fuselage mounted engines
        HYDPR    - Hydraulic system pressure                        [psi]  #default value equal to 3000
        VARSWP   - Wing variable sweep weight penalty factor               # 0 for fixed-geometry wing
                                                                           # 1 for full variable-sweep wing

    Outputs:
      WHYD      - Weight of hydraulic system group                  [kg]
    Properties Used: 
        N/A
    """ 
    # unpack inputs
    FPAREA = FPAREA/Units.ft**2 
    SW     = SW /Units.ft**2         # Convert kg to lb

    FPAREA = FPAREA * NFUSE
    
    #Calculate hydraulics weight using FLOPS methods
    WHYD = 0.57*(FPAREA+0.27*SW)*(1+0.03*FNEW+0.05*FNEF)*(3000/HYDPR)**0.35*(1+0.04*VARSWP)*VMAX**0.33
    
    WHYD = WHYD*Units.lb
    
    return WHYD
