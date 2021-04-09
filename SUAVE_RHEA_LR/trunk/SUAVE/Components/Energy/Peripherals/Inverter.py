## @ingroup Components-Energy-Peripherals
# Payload.py
# 
# Created:  Sep 2020, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE
from SUAVE.Components.Energy.Energy_Component import Energy_Component

# ----------------------------------------------------------------------
#  Payload Class
# ----------------------------------------------------------------------  
## @ingroup Components-Energy-Peripherals
class Inverter(Energy_Component):
    """A class representing a payload.
    
    Assumptions:
    None
    
    Source:
    N/A
    """          
    def __defaults__(self):
        """This sets the default efficiency and power-to weight ratios

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """            
        self.efficiency      = 0.9
        self.power_to_weight = 0.0 
        
