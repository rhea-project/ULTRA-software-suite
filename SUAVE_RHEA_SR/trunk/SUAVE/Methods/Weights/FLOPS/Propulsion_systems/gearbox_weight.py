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

def weight_gearbox(propeller,gearbox,motor,num_eng,power):
    
    """ Calculate the gearbox weight based on the NASA methods
    
   
    Source: 
        N/A
        
    Inputs:

    Outputs:
                                                    
      
    Assumptions:

    """ 
  

    # unpack inputs
    motor_RPM = motor.RPM
    prop_RPM  = propeller.RPM
    gear_year = gearbox.year

    if gear_year == '1980':
        Kgear = 43
    elif gear_year == '2000':
        Kgear = 34
    elif gear_year == 'future':
        Kgear = 26
    else:
        print('Warning! The gearbox year parameter is not found. A default value of 2000 is used')
        Kgear = 34

    WGB = Kgear * (power**0.76*motor_RPM**0.13)/(prop_RPM)**0.89 * num_eng   

    return WGB
      



