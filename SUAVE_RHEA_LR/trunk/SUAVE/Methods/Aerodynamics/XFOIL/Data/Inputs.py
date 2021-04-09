## @ingroup Methods-Aerodynamics-Xfoil-Data
# Inputs.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data

# ------------------------------------------------------------
#   Configuration
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Xfoil-Data
class Inputs(Data):
	""" A data class defining filenames for the Xfoil executable

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
		""" Defines the data structure  and defaults of airfoil configurations and cases 
	
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
		self.configuration  = Data()
		self.aircraft       = Data()
		self.cases          = Data()
		self.avl_bin_path   = r'C:\AE Software\Xfoil\xfoil.exe'
		
		filenames           = Data()
		filenames.geometry  = 'airfoil.dat'
		filenames.results   = []
		filenames.cases     = 'airfoil.cases'
		filenames.deck      = 'Xfoil_commands.run'
		filenames.reference_path = SUAVE.__path__[0] + '/temporary_files/'
		self.input_files = filenames
