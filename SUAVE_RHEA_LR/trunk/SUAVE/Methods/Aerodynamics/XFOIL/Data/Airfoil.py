## @ingroup Methods-Aerodynamics-Xfoil-Data
#Wing.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 
#  

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data

# ------------------------------------------------------------
#   Wing
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-AVL-Data
class Airfoil(Data):
	""" A class that defines parameters of the Xfoil airfoil

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

		self.tag = 'airfoil'

		self.sections           = Data()
		self.configuration      = Data()
		self.control_surfaces   = Data()

		self.configuration.nchordwise   = 240
		self.configuration.cspace       = 1.0


	def append_section(self,section):                       # CHANGE IF NEEDED!!!!
		""" adds a segment to the wing """

		# assert database type
		if not isinstance(section,Data):
			raise Exception('input component must be of type Data()')

		# store data
		self.sections.append(section)
		return


#### CHANGE IF NEEDED (FOR PLAN FLAPS OR CONTROL SURFACES

# ------------------------------------------------------------
#  AVL Wing Sections
# ------------------------------------------------------------

class Section(Data):
	""" A class that defines the sections of the aircraft wing in AVL.
	Each section can be thought of as a trapezoid

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
		""" Sets the defaunts of the aircraft wing geometry 
	
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
		self.tag    = 'section'
		self.origin = [0.0,0.0,0.0]
		self.chord  = 0.0
		self.twist  = 0.0
		self.airfoil_coord_file = None
		self.control_surfaces = Data()
		
				
	def append_control_surface(self,control):
		""" Adds a control_surface to the wing section in AVL
	
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

		# assert database type
		if not isinstance(control,Data):
			raise Exception('input component must be of type Data()')

		# store data
		self.control_surfaces.append(control)
		return


# ------------------------------------------------------------
#  AVL Control Surface
# ------------------------------------------------------------

class Control_Surface(Data):
	""" A class that defines the control surface geometry and deflection
	on the aircraft wing in AVL

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
		""" Sets the defaults of the control surface on the aircraft wing
		in AVL
	
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
		self.tag            = 'control_surface'
		self.gain           = 0.0
		self.x_hinge        = 0.0
		self.hinge_vector   = '0. 0. 0.'
		self.sign_duplicate = '+1'	# sign_duplicate: 1.0 or -1.0 - the sign of
						# the duplicate control on the mirror wing.
						# Use 1.0 for a mirrored control surface,
						# like an elevator. Use -1.0 for an aileron.

