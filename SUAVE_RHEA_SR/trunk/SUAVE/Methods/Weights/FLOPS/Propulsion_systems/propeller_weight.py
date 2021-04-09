## NASA Flops weight methods

#
# Created:  Mar 2021, S. Karpuk
# Modified: 
#           

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import numpy  as np
from SUAVE.Core import Units
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#  Engine Weight Calculation
# ----------------------------------------------------------------------

def weight_propeller(propeller,power,M,Neng):
    
    """ Calculate the propeller weight  based on the NASA methods
    
   
    Source: 
        N/A
        
    Inputs:

    Outputs:
                                                    
      
    Assumptions:
        The propeller activity factor is calculated using the propeller MAC
        double-acting propellers

    """ 
    # unpack inputs
    R     = propeller.tip_radius 
    R0    = propeller.hub_radius
    sigma = propeller.solidity
    Nb    = propeller.number_blades 
    Mat   = propeller.material 
    
    # Convert all values to the US units
    D    = 2*R / Units.feet
    cbar = (sigma*np.pi*R/Nb) / Units.feet
    SHP  = power*1.34102

    # Calculate the propeller weight
    if Mat == 'aluminium':
        Kw = 335
    elif Mat == 'composite':
        Kw = 180
    else:
        print('The propeller material is not specified! Assumed an aluminium')
    
    AF = 1562.5*cbar/D*(1-(R0/R)**4)

    WPROP = Kw * ((D/10)**2*(Nb/4)**0.7*(AF/100)**0.75*(Nb*D/20000)**0.5*(M+1)**0.5*(SHP/10*D**2)**0.12) * Neng
      
    WPROP = WPROP * Units.lb

    return WPROP
      



