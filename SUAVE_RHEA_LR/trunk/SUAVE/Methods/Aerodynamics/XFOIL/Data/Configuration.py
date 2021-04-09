## @ingroup Methods-Aerodynamics-Xfoil-Data
# Configuration.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data


# ------------------------------------------------------------
#   Configuration
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-AVL-Data
class Configuration(Data):
	"""A data class defining the reference parameters of the airfoil geoemtry and 
	flight configuration 

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
		""" Defines the data structure and default properties of the airfoil
		in Xfoil

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
		self.tag = 'configuration'
		self.parasite_drag                = 0.0
		
		self.reference_values             = Data()
		self.reference_values.cref        = 0.0
		self.reference_values.cg_coords   = [0.,0.,0.]

