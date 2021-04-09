## @ingroup Methods-Aerodynamics-Xfoil-Data
#Settings.py
# 
# Created:  Jul, S. Karpuk
# Modified: 
# 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Data
from .Cases import Run_Case

# ------------------------------------------------------------
#   Configuration
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Xfoil-Data
class Settings(Data):
        """ A class that defines important settings that call the Xfoil executable in addition to the 
        format of the result, batch and geometry filenames
        
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
                """ Defines naming convention for files created/used by AVL to compute analysus
        
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
                self.run_cases                = Run_Case.Container()
                self.filenames                = Data()
                self.discretization           = Data()
                self.number_panels            = 0
                
                self.discretization.defaults  = Data()
                self.discretization.defaults  = Xfoil_Discretization_Settings()

                self.filenames.xfoil_bin_name  = r'C:\AE_Software\Xfoil\xfoil.exe' # to call avl from command line. If avl is not on the system path, include absolute path to the avl binary
                self.filenames.run_folder      = 'xfoil_files' # local reference, will be attached to working directory from which avl was created
                self.filenames.features        = 'airfoil.xfoil'
                self.filenames.batch_template  = 'batch_{0:03d}.run'
                self.filenames.deck_template   = 'commands_{0:03d}.deck'
                self.filenames.output_template = 'results_{}.txt'
                self.filenames.case_template   = 'case_{0:03d}_{1:02d}'
                self.filenames.log_filename    = 'xfoil_log.txt'
                self.filenames.err_filename    = 'xfoil_err.txt'


# ------------------------------------------------------------
#  Xfoil Case
# ------------------------------------------------------------
## @ingroup Components-Airfoil
class Xfoil_Discretization_Settings(Data):
        """ A class that defines discretization of vortices on the airfoil
        
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
                """ Defines the spacing of vortices on the aifoil in Xfoil
               
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
                self.chordwise_vortices = 240                      

                
