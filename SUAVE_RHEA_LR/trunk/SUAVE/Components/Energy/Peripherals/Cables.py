## @ingroup Components-Energy-Peripherals
# Payload.py
# 
# Created:  Oct 2020, S. Karpuk
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
class Cables(Energy_Component):
    """A class representing cables.
    
    Assumptions:
    None
    
    Source:
    N/A
    """          
    def __defaults__(self):
        """This sets the default cable density.

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
        self.density = 0.0   # kg/m
        

