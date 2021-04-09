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
class TMS(Energy_Component):
    """A class representing a thermal managment system.
    
    Assumptions:
    None
    
    Source:
    N/A
    """          
    def __defaults__(self):
        """This sets the default TMS power density.

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
        self.power_density = 0.0   # kW/kg
        

