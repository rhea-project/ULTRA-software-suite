## @ingroup Analyses-Aerodynamics
# SU2_inviscid.py
#
# Created:  Sep 2016, E. Botero
# Modified: Jan 2017, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE imports
import SUAVE
from SUAVE.Core import Data, Units

# Local imports
from .Aerodynamics import Aerodynamics
from SUAVE.Input_Output.SU2.call_SU2_CFD import call_SU2_CFD
from SUAVE.Input_Output.SU2.write_SU2_cfg import write_SU2_cfg
from sklearn.gaussian_process.kernels import ExpSineSquared

# Package imports
import numpy as np
from contextlib import redirect_stdout
import io
import time
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
class SU2_inviscid(Aerodynamics):
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
        """Evaluates lift and drag using available surrogates.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        state.conditions.
          mach_number      [-]
          angle_of_attack  [radians]

        Outputs:
        inviscid_lift      [-] CL
        inviscid_drag      [-] CD

        Properties Used:
        self.surrogates.
          lift_coefficient [-] CL
          drag_coefficient [-] CD
        """  
        # Unpack
        surrogates = self.surrogates        
        conditions = state.conditions
        
        mach = conditions.freestream.mach_number
        AoA  = conditions.aerodynamics.angle_of_attack
        lift_model = surrogates.lift_coefficient
        drag_model = surrogates.drag_coefficient

        x = np.hstack((AoA,mach))    
        
        # Inviscid lift
        data_len = len(AoA)
        inviscid_lift = np.zeros([data_len,1])
        #for ii,_ in enumerate(AoA):
            #inviscid_lift[ii] = lift_model.predict([np.array([AoA[ii][0],mach[ii][0]])]) #sklearn fix
            #inviscid_lift[ii] = lift_model.predict_values([np.array([AoA[ii][0],mach[ii][0]])])[:, 0]

        trap = io.StringIO()
        with redirect_stdout(trap):
            inviscid_lift = np.transpose([lift_model.predict_values(x)[:, 0]])
            
        conditions.aerodynamics.lift_breakdown.inviscid_wings_lift       = Data()
        conditions.aerodynamics.lift_breakdown.inviscid_wings_lift.total = inviscid_lift
        state.conditions.aerodynamics.lift_coefficient                   = inviscid_lift
        state.conditions.aerodynamics.lift_breakdown.compressible_wings  = inviscid_lift

        # Inviscid drag, zeros are a placeholder for possible future implementation
        inviscid_drag                                              = np.zeros([data_len,1])       
        state.conditions.aerodynamics.inviscid_drag_coefficient    = inviscid_drag
        
        return inviscid_lift, inviscid_drag


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
          coefficients     [-] CL and CD
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
        
        AoA  = training.angle_of_attack
        mach = training.Mach 
        CL   = np.zeros([len(AoA)*len(mach),1])
        CD   = np.zeros([len(AoA)*len(mach),1])

        # Condition input, local, do not keep (k is used to avoid confusion)
        konditions              = Data()
        konditions.aerodynamics = Data()

        if self.training_file is None:
            # Calculate aerodynamics for table
            table_size = len(AoA)*len(mach)
            xy = np.zeros([table_size,2])
            count = 0
            time0 = time.time()
            for i,_ in enumerate(AoA):
                for j,_ in enumerate(mach):
                    
                    xy[count,:] = np.array([AoA[i],mach[j]])
                    # Set training conditions
                    konditions.aerodynamics.angle_of_attack = AoA[i]
                    konditions.aerodynamics.mach            = mach[j]
                    
                    CL[count],CD[count] = call_SU2(konditions, settings, geometry)
                    count += 1
            
            time1 = time.time()
            
            print('The total elapsed time to run SU2: '+ str(time1-time0) + '  Seconds')
        else:
            data_array = np.loadtxt(self.training_file)
            xy         = data_array[:,0:2]
            CL         = data_array[:,2:3]
            CD         = data_array[:,3:4]

        # Save the data
        np.savetxt(geometry.tag+'_data.txt',np.hstack([xy,CL,CD]),fmt='%10.8f',header='AoA Mach CL CD')

        # Store training data
        training.coefficients = np.hstack([CL,CD])
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
          coefficients     [-] CL and CD
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
        CL_data   = training.coefficients[:,0]
        CD_data   = training.coefficients[:,1]
        xy        = training.grid_points

        #xlimits   = np.array([[np.amin(AoA_data) - 0.01745,np.amax(AoA_data) + 0.01745], \
        #                     [np.amin(CL_data) - 0.1, np.amax(CL_data) + 0.5]])      
        xlimits   = np.array([[np.amin(xy[:,0]) - 0.1,np.amax(xy[:,0]) + 0.01745],[-1, np.amax(CL_data) + 0.5]])   

        # RMTB Method
        cl_surrogate = RMTB(num_ctrl_pts=60, xlimits=xlimits, nonlinear_maxiter=100, energy_weight=1e-12)
        cd_surrogate = RMTB(num_ctrl_pts=60, xlimits=xlimits, nonlinear_maxiter=100, energy_weight=1e-12)
        cl_surrogate.set_training_values(xy, np.transpose([CL_data]))
        cd_surrogate.set_training_values(xy, np.transpose([CD_data]))
        cl_surrogate.train()
        cd_surrogate.train()
       
        # Gaussian Process New
        #gp_kernel_ES = ExpSineSquared(length_scale=1.0, periodicity=1.0, length_scale_bounds=(1e-5,1e5), periodicity_bounds=(1e-5,1e5))
        #regr_cl = gaussian_process.GaussianProcessRegressor(kernel=gp_kernel_ES)
        #regr_cd = gaussian_process.GaussianProcessRegressor(kernel=gp_kernel_ES)
        #cl_surrogate = regr_cl.fit(xy, CL_data)
        #cd_surrogate = regr_cd.fit(xy, CD_data)  
        
        # KNN
        #regr_cl = neighbors.KNeighborsRegressor(n_neighbors=1,weights='distance')
        #regr_cd = neighbors.KNeighborsRegressor(n_neighbors=1,weights='distance')
        #cl_surrogate = regr_cl.fit(xy, CL_data)
        #cd_surrogate = regr_cd.fit(xy, CD_data)  
        
        # SVR
        #regr_cl = svm.SVR(C=500.)
        #regr_cd = svm.SVR()
        #cl_surrogate = regr_cl.fit(xy, CL_data)
        #cd_surrogate = regr_cd.fit(xy, CD_data)          
        
        
        self.surrogates.lift_coefficient = cl_surrogate
        self.surrogates.drag_coefficient = cd_surrogate
        
        # Standard subsonic test case
        AoA_points = np.linspace(-7.,7.,100)*Units.deg
        mach_points = np.linspace(.27,.9,100)      

        num_a = len(AoA_points)
        num_M = len(mach_points)
        
        AoA_mesh,mach_mesh = np.meshgrid(AoA_points,mach_points)
        
        CL_sur = np.zeros(np.shape(AoA_mesh))
        CD_sur = np.zeros(np.shape(AoA_mesh))        

        x = np.zeros((num_a,num_M, 2))
        x[:, :, 0] = np.outer(np.linspace(-8.0*3.14/180, 8.0*3.14/180, num_a), np.ones(num_M))
        x[:, :, 1] = np.outer(np.ones(num_a), np.linspace(0.29, 0.86, num_M))
        
        CL_sur = cl_surrogate.predict_values(x.reshape((num_a * num_M, 2)))[:, 0].reshape((num_a, num_M))
        CD_sur = cd_surrogate.predict_values(x.reshape((num_a * num_M, 2)))[:, 0].reshape((num_a, num_M))

        '''for jj in range(len(AoA_points)):
            for ii in range(len(mach_points)):
                CL_sur[ii,jj] = cl_surrogate.predict([np.array([AoA_mesh[ii,jj],mach_mesh[ii,jj]])])
                CD_sur[ii,jj] = cd_surrogate.predict([np.array([AoA_mesh[ii,jj],mach_mesh[ii,jj]])])  #sklearn fix       

        fig = plt.figure('Coefficient of Lift Surrogate Plot')    
        plt_handle = plt.contourf(AoA_mesh/Units.deg,mach_mesh,CL_sur,levels=None)
        #plt.clabel(plt_handle, inline=1, fontsize=10)
        cbar = plt.colorbar()
        plt.scatter(xy[:,0]/Units.deg,xy[:,1])
        plt.xlabel('Angle of Attack (deg)')
        plt.ylabel('Mach Number')
        cbar.ax.set_ylabel('Coefficient of Lift')'''

        plt.clf()
        plt.plot( xy[:,1], xy[:,0]*180/3.14, "o")
        plt.contour(x[:, :, 1], x[:, :, 0] * 180/3.14, CL_sur, 20)
        plt.pcolormesh(x[:, :, 1], x[:, :, 0] * 180/3.14, CL_sur, cmap=plt.get_cmap("rainbow"))
        plt.xlabel("Mach number")
        plt.ylabel("Angle of Attack (deg)")
        plt.title("CL surrogate")
        plt.colorbar()        

        plt.savefig('surrogatelift.png')

        return



# ----------------------------------------------------------------------
#  Helper Functions
# ----------------------------------------------------------------------

def call_SU2(conditions,settings,geometry):
    """Calculates lift and drag using SU2

    Assumptions:
    None

    Source:
    N/A

    Inputs:
    conditions.
      mach_number        [-]
      angle_of_attack    [radians]
    settings.
      half_mesh_flag     <boolean> Determines if a symmetry plane is used
      parallel           <boolean>
      processors         [-]
      maximum_iterations [-]
    geometry.
      tag
      reference_area     [m^2]

    Outputs:
    CL                   [-]
    CD                   [-]

    Properties Used:
    N/A
    """      

    half_mesh_flag = settings.half_mesh_flag
    tag            = geometry.tag
    parallel       = settings.parallel
    processors     = settings.processors 
    iters          = settings.maximum_iterations
    
    SU2_settings = Data()
    if half_mesh_flag == False:
        SU2_settings.reference_area  = geometry.reference_area
    else:
        SU2_settings.reference_area  = geometry.reference_area/2.
    SU2_settings.mach_number     = conditions.aerodynamics.mach
    SU2_settings.angle_of_attack = conditions.aerodynamics.angle_of_attack / Units.deg
    SU2_settings.maximum_iterations = iters
    
    # Build SU2 configuration file
    write_SU2_cfg(tag, SU2_settings)
    
    # Run SU2
    CL, CD = call_SU2_CFD(tag,parallel,processors)
        
    return CL, CD
