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
#  unusable fuel Weight Calculation
# ----------------------------------------------------------------------

def weight_unusable_fuel(ac_type,num_eng,rat_thrust,wing_area,num_tank,ac_max_fuel_cap):
    
    """ Calculate the unusable fuel weight  based on the FLOPS methods    
   
    Source: 
        N/A
        
    Inputs:
        ac_type - aircraft type                                             [Unitless]
        num_eng - number of engines                                         [Unitless]
        
        rat_thrust - Rated thrust of each scaled engine                     [N]
        wing_area  - Reference wing area                                    [m²]
        ac_max_fuel_cap - Aircraft maximum fuel capacity                    [kg]
        num_tank - number of fuel tanks
ns
    Outputs:
      WUF - Weight of unusable fuel                                         [kg]                                             
      
    Properties Used:
        N/A
    """ 
     # unpack inputs
    rat_thrust      = rat_thrust*0.22481 #N to lb
    wing_area       = wing_area /Units.ft**2  
    ac_max_fuel_cap = ac_max_fuel_cap/Units.lb  # Convert meter to ft

    #Calculate the unusable fuel  weight using FLOPS methods
    if ac_type=='transport':
        WUF = 11.5*num_eng*rat_thrust**0.2+0.007*wing_area+1.6*num_tank*ac_max_fuel_cap**0.28
    else:
        WUF = 11.5*num_eng*rat_thrust**0.2+0.04*wing_area

    WUF = WUF*Units.kg 
    return WUF

