## @ingroup methods-Aerodynamics-Quasi_3D-Drag-Q3D
# AVL.py
#
# Created: Jul 2019, S. Karpuk 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data
from SUAVE.Core import redirect

from SUAVE.Analyses import Process
from SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics                   import Aerodynamics
from SUAVE.Analyses.Mission.Segments.Conditions.Conditions                     import Conditions
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Drag                      import  parasite_drag_wing,compressibility_drag_wing
from SUAVE.Methods.Aerodynamics.Quasi_3D.Quasi_3D_wing                         import Quasi_3D_wing
from SUAVE.Methods.Geometry.Three_Dimensional.compute_wing_section_coordinates import compute_wing_section_coordinates as wing_sections
from SUAVE.Methods.Aerodynamics.XFOIL.Data.Settings                            import Settings
from SUAVE.Methods.Aerodynamics.Quasi_3D.Data.Settings                         import Settings as freestream_Settings

# Package imports
import numpy as np
import os
from sklearn import gaussian_process

class Q3D(Aerodynamics):
    """This builds a surrogate and computes drag using AVL and Xfoil.

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
        self.tag                             = 'Q3D'
        self.keep_files                      = True            
           
        self.settings                        = Settings()
        self.fstr_settings                   = freestream_Settings()

        # Define Q3D paramters
        self.number_of_drag_setions          = self.fstr_settings.discretization.defaults.number_of_Xfoil_drag_setions
        self.ref_line_percent                = self.fstr_settings.discretization.defaults.ref_line_percent  
        
        # Conditions table, used for surrogate model training
        self.training                        = Data()
        
        # Standard subsonic/transonic aircarft
        self.training.angle_of_attack         = self.fstr_settings.training.angle_of_attack
        self.training.Mach                    = self.fstr_settings.training.Mach   
        self.training.Re                      = self.fstr_settings.training.Re

        # Surrogate model
        self.surrogates                      = Data()
        
        # Regression Status
        self.regression_flag                 = False

        
    def initialize(self):
        """Drives functions to get training samples and build a drag surrogate.

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
        self.tag     = 'Q3D_analysis_of_{}'.format(geometry.tag)
        run_folder   = self.settings.filenames.run_folder

        # Sample training data
        self.sample_training()

        # Build surrogate
        self.build_surrogate()
        
        return
        
    def evaluate(self,state,settings,geometry):
        """Computes the parasite drag due to wings using a Quasi-3D approach for the Main wing
            and the Fidelity Zero for any other wing-like component

        Assumptions:
        Basic fit

        Source:

        Inputs:
        state.conditions.freestream.
          mach_number                                [Unitless]
          temperature                                [K]
          reynolds_number                            [Unitless]
        geometry.
          areas.reference                            [m^2]
          chords.mean_aerodynamic                    [m]
          thickness_to_chord                         [Unitless]
          sweeps.quarter_chord                       [radians]
          aspect_ratio                               [Unitless]
          spans.projected                            [m]
          areas.exposed                              [m^2]
          areas.affected                             [m^2]
          areas.wetted                               [m^2]
          transition_x_upper                         [Unitless]
          transition_x_lower                         [Unitless]
          
          
        Outputs:
        wing_parasite_drag                           [Unitless]

        Properties Used:
        N/A
        
        """
        # Unpack
        surrogates    = self.surrogates       
        conditions    = state.conditions

        mach          = conditions.freestream.mach_number
        AoA           = conditions.aerodynamics.angle_of_attack
        Re            = conditions.freestream.reynolds_number
        drag_model    = surrogates.drag_coefficient

        data_len       = len(AoA)
        viscous_drag   = np.zeros([data_len,1])
        cd_c           = np.zeros([data_len,1])
        
        # Run analysis depending on the aircraft component
        if geometry.tag == 'main_wing':

            # Run Q3D
            for ii,_ in enumerate(AoA):
                viscous_drag[ii] = drag_model.predict([np.array([AoA[ii][0],mach[ii][0],Re[ii][0]/10E7])])

            # dump data to conditions
            wing_result             = Data(parasite_drag_coefficient = viscous_drag)
            wing_compressibe_result = Data(compressibility_drag      = cd_c)

            # Store viscous and compressible drag results
            state.conditions.aerodynamics.drag_breakdown.parasite[geometry.tag]     = wing_result
            state.conditions.aerodynamics.drag_breakdown.compressible[geometry.tag] = wing_compressibe_result
            '''f = open("results.dat","w")
            f.write(str(wing_result))
            f.close()'''
        else:
            # Run Fidelity zero
            parasite_drag_wing(state,settings,geometry)
            compressibility_drag_wing(state,settings,geometry)

        return 

    def sample_training(self):
        """Call methods to run Q3D for sample point evaluation.

        Assumptions:
        Returned drag values are not meaningful.

        Source:
        N/A

        Inputs:
        see properties used

        Outputs:
        self.training.
          coefficients     [-] CD
          grid_points      [radians,-] angles of attack and mach numbers 

        Properties Used:
        self.training.
          Re               [-]
          angle_of_attack  [radians]
          Mach             [-]
        """
        
        # Unpack
        training = self.training

        AoA  = self.training.angle_of_attack        
        Mach = self.training.Mach                   
        Re   = self.training.Re

        #Remin = 1.225*343*0.1/(1.81E-5)/(6)
        #Remax = 1.225*343*0.55/(1.81E-5)/(0.5*6)

        # Calculate aerodynamics for table
        table_size = len(AoA)*len(Mach)*len(Re)
        xy         = np.zeros([table_size,3])  

        for i,_ in enumerate(Mach):
            for j,_ in enumerate(AoA):
                for k,_ in enumerate(Re):
                    xy[i*len(Re)*len(AoA)+j*len(Re)+k,:] = np.array([AoA[j],Mach[i],Re[k]/10E7])

        # Run Q3D to obdatin the values
        wing_parasite_drag = Quasi_3D_wing(self)

        # Store training data
        training.coefficients = wing_parasite_drag
        training.grid_points  = xy

        return


    def build_surrogate(self):
        """Builds a surrogate based on sample evalations using a Guassian process.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        self.training.
          coefficients             [-] CD 
          grid_points              [radians,-] angles of attack,mach and Reynolds numbers 

        Outputs:
        self.surrogates.
          drag_coefficient       <Guassian process surrogate>

        Properties Used:
        No others
        """
        
        # Unpack data
        training                         = self.training
        AoA_data                         = training.angle_of_attack
        mach_data                        = training.Mach
        Re_data                          = training.Re
        CD_data                          = training.coefficients
        xy                               = training.grid_points

        # Gaussian Process New
        regr_cd                          = gaussian_process.GaussianProcessRegressor()
        cd_surrogate                     = regr_cd.fit(xy, CD_data)
        self.surrogates.drag_coefficient = cd_surrogate

        return
