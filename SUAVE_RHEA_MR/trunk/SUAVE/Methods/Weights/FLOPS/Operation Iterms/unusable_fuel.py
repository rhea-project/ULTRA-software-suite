## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

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

def weight_unusable_fuel(ac_type,
                         num_eng,
                         rat_thrust,
                         wing_area,
                         num_tank,
                         ac_max_fuel_cap):
    """ Calculate the unusable fuel weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type - aircraft type                 : 1 for transport/general aviation aircraft
                                                  2 for fighter/attack aircraft
        num_eng - number of engines
        
        rat_thrust - Rated thrust of each scaled engine                     [kN]
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
    rat_thrust=rat_thrust*224.81 #kN to lb
    wing_area=wing_area /Units.ft**2  
    ac_max_fuel_cap=ac_max_fuel_cap/Units.lb  # Convert meter to ft

    #Calculate the unusable fuel  weight using FLOPS methods
    WUF=11.5*num_eng*rat_thrust**0.2+0.007*wing_area+1.6*num_tank*ac_max_fuel_cap**0.28 if ac_type==1  else 11.5*num_eng*rat_thrust**0.2+0.04*wing_area
    WUF=WUF*Units.kg 
    return WUF

