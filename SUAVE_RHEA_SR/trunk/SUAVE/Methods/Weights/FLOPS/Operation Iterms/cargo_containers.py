## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


import pandas as pd
import math as m
from SUAVE.Core import Units

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Cargo containers Weight Calculation
# ----------------------------------------------------------------------

def weight_cargo_containers(weight_cargo):
    """ Calculate the weight of the cargo containers  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
      weight_cargo  - Weight of the cargo that will be placed in containers  [kg]
       
ns
    Outputs:
     weight_cargo_containers - Weight of cargo containers                  [kg]
    Properties Used:
        N/A
    """ 
 
    weight_cargo=weight_cargo/Units.lb #convert kg to lb
    num_cargo=m.ceil(weight_cargo/950) #assuming that each cargo container contains 9520lb of cargo
    weight_cargo_containers=175*num_cargo #Each cargo container is assumed to weigh 157lb
    weight_cargo_containers=weight_cargo_containers*Units.kg
    return weight_cargo_containers
    
