## @ingroup Analyses-Aerodynamics
# AVL_Inviscid_loads.py
#
# Created:  Apr 2017, M. Clarke 
# Modified: Jan 2018, W. Maier
#           Oct 2018, M. Clarke
#           Jan 2020, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports
import SUAVE
from SUAVE.Core import Units, Data
from SUAVE.Core import redirect

from SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics import Aerodynamics
from SUAVE.Analyses.Mission.Segments.Conditions.Conditions   import Conditions

from SUAVE.Methods.Aerodynamics.AVL.write_geometry         import write_geometry
from SUAVE.Methods.Aerodynamics.AVL_Loads.write_run_cases  import write_run_cases
from SUAVE.Methods.Aerodynamics.AVL_Loads.write_input_deck import write_input_deck
from SUAVE.Methods.Aerodynamics.AVL_Loads.run_analysis     import run_analysis
from SUAVE.Methods.Aerodynamics.AVL_Loads.translate_data   import translate_conditions_to_cases, translate_results_to_conditions
from SUAVE.Methods.Aerodynamics.AVL.purge_files            import purge_files
from SUAVE.Methods.Aerodynamics.AVL_Loads.Data.Settings    import Settings
from SUAVE.Methods.Aerodynamics.AVL.Data.Cases             import Run_Case

# Package imports
import time
import pylab as plt
import os
import sklearn
from sklearn import gaussian_process
import numpy as np
import sys
from shutil import rmtree
from warnings import warn

# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class AVL_Inviscid_loads(Aerodynamics):
    """This script obtains loads from the AVL solution.

    Assumptions:
    None

    Source:
    None
    """     
    def __defaults__(self):
        """This sets the default values and methods for the analysis.

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
        self.tag                             = 'avl'
        self.keep_files                      = True
                
        self.current_status                  = Data()        
        self.current_status.batch_index      = 0
        self.current_status.batch_file       = None
        self.current_status.deck_file        = None
        self.current_status.cases            = None      
        self.geometry                        = None   
        
        self.settings                        = Settings()
        self.settings.filenames.log_filename = sys.stdout
        self.settings.filenames.err_filename = sys.stderr        
        self.settings.spanwise_vortices      = None 
        self.settings.chordwise_vortices     = None         
        
        # Conditions table, used for surrogate model training
        self.training                        = Data()
        
        # Regression Status
        self.regression_flag                 = False

    def initialize(self,state,spanwise_vortices,chordwise_vortices):
        """Drives functions evaluate the AVL-case to determine loads.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        self.tag = 'avl_analysis_of_{}'.format(geometry.tag)

        Properties Used:
        self.geometry.tag
        """  
        geometry     = self.geometry
        self.tag     = 'avl_analysis_of_{}'.format(geometry.tag)
        run_folder   = self.settings.filenames.run_folder
        
        # check if user specifies number of spanwise vortices
        if spanwise_vortices == None:
            pass
        else:
            self.settings.discretization.defaults.wing.spanwise_vortices = spanwise_vortices  
        
        # check if user specifies number of chordise vortices 
        if chordwise_vortices == None:
            pass
        else:
            self.settings.discretization.defaults.wing.chordwise_vortices = chordwise_vortices

        # obtain loads from AVL
        self.loads_analysis(state)
    
        return
        

    def loads_analysis(self,state):
        """Call methods to run AVL for loads analysis.

        Assumptions:
        Returned drag values are not meaningful.

        Source:
        N/A

        Inputs:
        see properties used

        Outputs:

        Properties Used:
        self.geometry.tag  <string>
        state.conditions.     
          aerodynamics.lift_coefficient [-]
          freestream.mach               [-]
        """
        
        # Unpack 
        geometry = self.geometry  
        
        CL      = state.conditions.aerodynamics.lift_coefficient
        mach    = state.conditions.freestream.mach  
        
        # Calculate aerodynamics for table 
        time0      = time.time()

        run_conditions = Aerodynamics()
        run_conditions.freestream.density           = 0     # Density not used in inviscid computation therefore set to zero. Used for dynamic analysis which is under development
        run_conditions.freestream.gravity           = 9.81        
        run_conditions.aerodynamics.CL              = CL
        run_conditions.freestream.mach_number       = mach
        run_conditions.aerodynamics.angle_of_attack = [0]
            
        #Run Analysis at AoA[i] and mach[j]
        results =  self.evaluate_conditions(run_conditions)
                              
        time1 = time.time()
        
        print('The total elapsed time to run AVL: '+ str(time1-time0) + '  Seconds')

        # Obtain shear and bending moments
        loads = self.read_loads()
                    
        # Store loads data
        self.loads = loads

        return        
           

# ----------------------------------------------------------------------
#  Helper Functions
# ----------------------------------------------------------------------
        
    def evaluate_conditions(self,run_conditions):
        """Process vehicle to setup geometry, condititon, and configuration.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        run_conditions <SUAVE data type> aerodynamic conditions; until input
                method is finalized, will assume mass_properties are always as 
                defined in self.features

        Outputs:
        results        <SUAVE data type>

        Properties Used:
        self.settings.filenames.
          run_folder
          output_template
          batch_template
          deck_template
        self.current_status.
          batch_index
          batch_file
          deck_file
          cases
        """           
        
        # unpack
        run_folder                       = os.path.abspath(self.settings.filenames.run_folder)
        loads_output_template            = self.settings.filenames.loads_output_template
        output_template                  = self.settings.filenames.output_template
        batch_template                   = self.settings.filenames.batch_template
        deck_template                    = self.settings.filenames.deck_template
        
        # rename defaul avl aircraft tag
        self.settings.filenames.features = self.geometry._base.tag + '.avl'
        
        # update current status
        self.current_status.batch_index += 1
        batch_index                      = self.current_status.batch_index
        self.current_status.batch_file   = batch_template.format(batch_index)
        self.current_status.deck_file    = deck_template.format(batch_index)
               
        # control surfaces
        num_cs = 0       
        for wing in self.geometry.wings:
            for segment in wing.Segments:
                wing_segment =  wing.Segments[segment]
                section_cs = len(wing_segment.control_surfaces)
                if section_cs != 0:
                    cs_shift = True
                num_cs =  num_cs + section_cs

        # translate conditions
        cases                            = translate_conditions_to_cases(self,run_conditions)        
        for case in cases:
            cases[case].stability_and_control.number_control_surfaces = num_cs

        self.current_status.cases        = cases 
        
        # case filenames
        for case in cases:
            cases[case].result_filename       = output_template.format(case)
            cases[case].loads_result_filename = loads_output_template.format(case)
    
        # write the input files
        with redirect.folder(run_folder,force=False):
            write_geometry(self)
            write_run_cases(self)
            write_input_deck(self)
    
            # RUN AVL!
            results_avl = run_analysis(self)
    
        # translate results
        results = translate_results_to_conditions(cases,results_avl)
    
        if not self.keep_files:
            rmtree( run_folder )
    
        return results


    def read_loads(self):
        """reads loads file into an array

        Assumptions:
        None

        Source:
        N/A

        Inputs:

        Outputs:
        loads        <SUAVE data type>
        
        """
        run_folder = os.path.abspath(self.settings.filenames.run_folder)
        N_strips   = self.settings.discretization.defaults.wing.spanwise_vortices + 2
        loads      = np.zeros((N_strips, 3))
        with redirect.folder(run_folder,force=False):
            with open('loads_results.txt') as f:
                for line in f:
                    if 'Surface:   1' in line:
                        for _ in range(4):
                            next(f)
                        for i in range(N_strips):
                            num_line = f.readline()
                            string = np.array(num_line.split())
                            loads[i,:] = string.astype(np.float)

        return loads
