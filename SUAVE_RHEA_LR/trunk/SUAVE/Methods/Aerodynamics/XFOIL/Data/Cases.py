## @ingroup Methods-Aerodynamics-Xfoil-Data
# Cases.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 
#  

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data
from SUAVE.Core import DataOrdered 

# ------------------------------------------------------------
#  Xfoil Case
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-AVL-Data
class Run_Case(Data):
    """ A data class defining the parameters for the analysis cases 
    including angle of attack and mach number 

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        None

    Outputs:
        None

    Properties Used:
        N/A
    """    
    
    def __defaults__(self):
        """Defines the data structure and defaults of aerodynamics coefficients 

        Assumptions:
            None
    
        Source:
            None
    
        Inputs:
            None
    
        Outputs:
            None
    
        Properties Used:
            N/A
        """ 

        self.index                      = 0		# Will be overwritten when passed to an AVL_Callable object
        self.tag                        = 'case'

        self.conditions                 = Data()
        self.flap                       = Data()
        free                            = Data()
        aero                            = Data()

        free.mach                       = 0.0
        free.velocity                   = 0.0
        free.density                    = 1.225
        free.gravitational_acceleration = 9.81

        aero.parasite_drag              = 0.0
        aero.angle_of_attack            = 0.0

        self.flap.deflections           = None
        self.conditions.freestream      = free
        self.conditions.aerodynamics    = aero

        self.result_filename            = None
        self.eigen_result_filename      = None
 

    def append_flap_deflection(self,flap_tag,deflection):
        """ Adds a control deflection case 

	Assumptions:
	    None
    
	Source:
	    None
    
	Inputs:
	    None
    
	Outputs:
	    None
    
	Properties Used:
	    N/A
	"""         
        flap_deflection              = Data()
        flap_deflection.tag          = flap_tag
        flap_deflection.deflection   = deflection
        if self.flap.flap_deflections is None:
            self.flap.flap_deflections = Data()
        self.flap.flap_deflections.append(flap_deflection)

        return

class Container(DataOrdered):
    """ A data class for the addition of a cases to the set of run cases

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        None

    Outputs:
        None

    Properties Used:
        N/A
    """    
    def append_case(self,case):
        """ Adds a case to the set of run cases "
        
	Assumptions:
	    None
    
	Source:
	    None
    
	Inputs:
	    None
    
	Outputs:
	    None
    
	Properties Used:
	    N/A
	"""         
        case.index = len(self)+1
        self.append(case)

        return
    
    
# ------------------------------------------------------------
#  Handle Linking
# ------------------------------------------------------------
Run_Case.Container = Container
