## @ingroup Analyses-Aerodynamics
# SU2_inviscid.py
#
# Created:  Nov 2019, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports
import SUAVE
from SUAVE.Core import Data, Units

# Local imports
from .Aerodynamics import Aerodynamics
from sklearn.gaussian_process.kernels import ExpSineSquared

# Package imports
import numpy as np
from contextlib import redirect_stdout
import io
import time
import glob
import pylab as plt
import sklearn
from sklearn import gaussian_process
from sklearn import neighbors
from sklearn import svm
from smt.surrogate_models import RMTB

# ----------------------------------------------------------------------
#  Class
# ----------------------------------------------------------------------
## @ingroup Analyses-Aerodynamics
class SU2_wave_drag(Aerodynamics):
    """This builds a surrogate and computes lift and drag using SU2

    Assumptions:
    Inviscid, subsonic

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
        self.tag = 'SU2_inviscid'

        self.geometry = Data()
        self.settings = Data()
        self.settings.half_mesh_flag     = True
        self.settings.parallel           = False
        self.settings.processors         = 1
        self.settings.maximum_iterations = 1500

        # Conditions table, used for surrogate model training
        self.training = Data()        
        self.training.angle_of_attack  = np.array([-2.,3.,8.]) * Units.deg
        self.training.Mach             = np.array([0.3,0.7,0.85])
        self.training.lift_coefficient = None
        self.training.drag_coefficient = None
        self.training_file             = None
        
        # Surrogate model
        self.surrogates = Data()
 
        
    def initialize(self):
        """Drives functions to get training samples and build a surrogate.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        None
        """
        trap = io.StringIO()
        with redirect_stdout(trap):

        # Sample training data
            self.sample_training()
                    
        # Build surrogate
            self.build_surrogate()


    def evaluate(self,state,settings,geometry):
        """Evaluates wave drag using available surrogates.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        state.conditions.
          mach_number      [-]
          angle_of_attack  [radians]

        Outputs:
        wave_drag      [-] CDwave

        Properties Used:
        self.surrogates.
          wave_drag_coefficient [-] CDwave
        """  
        # Unpack
        surrogates = self.surrogates        
        conditions = state.conditions
        
        mach = conditions.freestream.mach_number
        AoA  = conditions.aerodynamics.angle_of_attack
        drag_model = surrogates.drag_coefficient

        x = np.hstack((AoA,mach)) 
        
        # Wave drag, zeros are a placeholder for possible future implementation
        data_len = len(AoA)
        wave_drag = np.zeros([data_len,1])
        
        #for ii,_ in enumerate(AoA):
            # wave_drag[ii] = drag_model.predict([np.array([AoA[ii][0],mach[ii][0]])]) #sklearn fix
            # wave_drag[ii] = drag_model.predict_values([np.array([AoA[ii][0],mach[ii][0]])])[:, 0]
        trap = io.StringIO()
        with redirect_stdout(trap):
            wave_drag = np.transpose([drag_model.predict_values(x)[:, 0]])
        
            # correct for possible negative values
        #for ii,_ in enumerate(AoA):
        #    if wave_drag[ii] < 0:
        #        wave_drag[ii] = 0

        main_wing_tag = self.geometry.wings.main_wing.tag
        conditions.aerodynamics.drag_breakdown.compressible[main_wing_tag] = Data()
        conditions.aerodynamics.drag_breakdown.compressible[main_wing_tag].compressibility_drag = wave_drag
        
        return wave_drag


    def sample_training(self):
        """Call methods to run SU2 for sample point evaluation.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        see properties used

        Outputs:
        self.training.
          coefficients     [-]  CDwave
          grid_points      [radians,-] angles of attack and mach numbers 

        Properties Used:
        self.geometry.tag  <string>
        self.training.     
          angle_of_attack  [radians]
          Mach             [-]
        self.training_file (optional - file containing previous AVL data)
        """               
        # Unpack
        geometry = self.geometry
        settings = self.settings
        training = self.training
        
        AoA    = training.angle_of_attack
        mach   = training.Mach 
        CDwave = np.zeros([len(AoA)*len(mach),1])

        # Condition input, local, do not keep (k is used to avoid confusion)
        konditions              = Data()
        konditions.aerodynamics = Data()

        if self.training_file is None:           
            print('No wave drag file exist. Check the input files')
        else:
            # Read compressibility rag data from the text file
            data_array = np.loadtxt(self.training_file)
            xy         = data_array[:,0:2]
            CL         = data_array[:,2:3]
            CDwave     = data_array[:,3:4]


        # Store training data
        training.coefficients = np.hstack([CL,CDwave])
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
          coefficients     [-] CD
          grid_points      [radians,-] angles of attack and mach numbers 

        Outputs:
        self.surrogates.
          lift_coefficient <Guassian process surrogate>
          drag_coefficient <Guassian process surrogate>

        Properties Used:
        No others
        """  
        # Unpack data
        training  = self.training
        AoA_data  = training.angle_of_attack
        mach_data = training.Mach
        CD_data   = training.coefficients[:,1]
        xy        = training.grid_points

        #xlimits   = np.array([[np.amin(AoA_data) - 0.01745,np.amax(AoA_data) + 0.01745], \
        #                      [np.amin(CD_data) - 0.1, np.amax(mach_data) + 0.5]])
        xlimits   = np.array([[np.amin(xy[:,0]) - 0.1,np.amax(xy[:,0]) + 0.01745], \
                              [np.amin(CD_data) - 0.4, np.amax(mach_data) + 0.5]])  

        # RMTB Method
        cd_surrogate = RMTB(num_ctrl_pts=60, xlimits=xlimits, nonlinear_maxiter=100, energy_weight=1e-12)
        cd_surrogate.set_training_values(xy, np.transpose([CD_data]))
        cd_surrogate.train()
              
        # Gaussian Process New
        #gp_kernel_ES = ExpSineSquared(length_scale=1.0, periodicity=1.0, length_scale_bounds=(1e-5,1e5), periodicity_bounds=(1e-5,1e5))
        #regr_cd = gaussian_process.GaussianProcessRegressor(kernel=gp_kernel_ES)
        #cd_surrogate = regr_cd.fit(xy, CD_data)  
        
        # KNN
        #regr_cl = neighbors.KNeighborsRegressor(n_neighbors=1,weights='distance')
        #regr_cd = neighbors.KNeighborsRegressor(n_neighbors=1,weights='distance')
        #cd_surrogate = regr_cd.fit(xy, CD_data)  
        
        # SVR
        #regr_cd = svm.SVR()
        #cd_surrogate = regr_cd.fit(xy, CD_data)          

        
        self.surrogates.drag_coefficient = cd_surrogate
        
        # Standard subsonic test case
        AoA_points = np.linspace(-7.,7.,100)*Units.deg
        mach_points = np.linspace(.27,.9,100)      

        num_a = len(AoA_points)
        num_M = len(mach_points)
        
        AoA_mesh,mach_mesh = np.meshgrid(AoA_points,mach_points)
        
        CD_sur = np.zeros(np.shape(AoA_mesh))        

        x = np.zeros((num_a,num_M, 2))
        x[:, :, 0] = np.outer(np.linspace(-8.0*3.14/180, 8.0*3.14/180, num_a), np.ones(num_M))
        x[:, :, 1] = np.outer(np.ones(num_a), np.linspace(0.29, 0.86, num_M))
        
        CD_sur = cd_surrogate.predict_values(x.reshape((num_a * num_M, 2)))[:, 0].reshape((num_a, num_M))

        plt.clf()
        plt.plot( xy[:,1], xy[:,0]*180/3.14, "o")
        plt.contour(x[:, :, 1], x[:, :, 0] * 180/3.14, CD_sur, 20)
        plt.pcolormesh(x[:, :, 1], x[:, :, 0] * 180/3.14, CD_sur, cmap=plt.get_cmap("rainbow"))
        
        plt.xlabel("Mach number")
        plt.ylabel("Angle of Attack (deg)")
        plt.title("CD_compressible surrogate")
        plt.colorbar()        

        plt.savefig('surrogatedrag.png')

        return

