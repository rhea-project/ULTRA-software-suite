## @ingroup Methods-Aerodynamics-AVL-Data
#Settings.py
# 
# Created:  Dec 2014, T. Momose
# Modified: Jan 2016, E. Botero
#           Oct 2018, M. Clarke

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Data, Units
from .Cases import Run_Case

import numpy as np

# ------------------------------------------------------------
#   Configuration
# ------------------------------------------------------------

## @ingroup Methods-Aerodynamics-AVL-Data
class Settings(Data):
        """ A class that defines important settings that call the AVL executable in addition to the 
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
                self.flow_symmetry            = Data()
                self.discretization           = Data()
                self.training                 = Data()
                self.number_control_surfaces     = 0
                
                self.discretization.defaults  = Data()
                self.discretization.surfaces  = Data()
                self.discretization.defaults.number_of_Xfoil_drag_setions     = 10                                # Usually, around 10
                self.discretization.defaults.ref_line_percent                 = 0.25                              # 0.25 (subsonic) or 0.5 (transonic)
                self.discretization.defaults.wing                             = AVL_Discretization_Settings()
                self.discretization.defaults.fuselage                         = AVL_Discretization_Settings()
                self.discretization.defaults.fuselage.spanwise_spacing        = 3
                self.discretization.defaults.fuselage.spanwise_spacing_scheme = 'equal'
                self.discretization.defaults.fuselage.nose_interpolation      = 'parabolic'
                self.discretization.defaults.fuselage.tail_interpolation      = 'linear'

                self.filenames.avl_bin_name            = r'C:\AE Software\AVL\avl.exe' # to call avl from command line. If avl is not on the system path, include absolute path to the avl binary
                self.filenames.run_folder              = 'Q3D_files'                            # local reference, will be attached to working directory from which avl was created
                self.filenames.features                = 'aircraft.avl'
                self.filenames.mass_file               = 'aircraft.mass'
                self.filenames.batch_template          = 'batch_{0:03d}.run'
                self.filenames.deck_template           = 'commands_{0:03d}.deck'
                self.filenames.output_template         = 'results_{}.txt'
                self.filenames.Trefftz_output_template = 'Trefftz_results_{}.txt'
                self.filenames.case_template           = 'case_{0:03d}_{1:02d}'
                self.filenames.log_filename            = 'avl_log.txt'
                self.filenames.err_filename            = 'avl_err.txt'

                # Discretization settings for Q3D
                self
                
                # Definition of freestream conditions for AVL and Xfoil 
                self.training.angle_of_attack  = np.array([-2., 0., 2., 5., 7., 10.])*Units.degrees
                self.training.Mach             = np.array([0.05,0.25,0.35,0.5])
                self.training.Re               = np.array([1E6,2E6,4E6,6E6])
                
                #--------------------------------------------------------------------------
                #           SUAVE-AVL dynamic stability analysis under development
                # self.filenames.stability_output_template = 'eigen_vals_results_{}.txt'
                #
                #--------------------------------------------------------------------------

                self.flow_symmetry.xz_plane         = 0	# Symmetry across the xz-plane, y=0
                self.flow_symmetry.xy_parallel      = 0 # Symmetry across the z=z_symmetry_plane plane
                self.flow_symmetry.z_symmetry_plane = 0.0



# ------------------------------------------------------------
#  AVL Case
# ------------------------------------------------------------
## @ingroup Components-Wings
class AVL_Discretization_Settings(Data):
        """ A class that defines discretization of vortices on the aircraft wing
        
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
                """ Defines the spacing of vortices on lifting surface in AVL
                SPACING SCHEMES:
                	- 'cosine' : ||  |    |      |      |    |  || (bunched at both ends)
                	- '+sine'  : || |  |   |    |    |     |     | (bunched at start)
                	- 'equal'  : |   |   |   |   |   |   |   |   | (equally spaced)
                	- '-sine'  : |     |     |    |    |   |  | || (bunched at end)
               
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
                self.chordwise_vortices        = 5
                self.chordwise_spacing_scheme  = 'equal'
                self.spanwise_vortices         = 10
                self.spanwise_spacing_scheme   = 'cosine'

