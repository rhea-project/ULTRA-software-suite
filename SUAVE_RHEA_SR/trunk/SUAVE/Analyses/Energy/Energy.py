## @ingroup Analyses-Energy
# Energy.py
#
# Created:  
# Modified: Feb 2016, Andrew Wendorff

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import numpy as np
import pylab as plt

from SUAVE.Analyses import Analysis
from smt.surrogate_models import RMTB
from SUAVE.Core import Data , Units

# ----------------------------------------------------------------------
#  Analysis
# ----------------------------------------------------------------------

## @ingroup Analyses-Energy
class Energy(Analysis):
    """ SUAVE.Analyses.Energy.Energy()
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
        self.tag     = 'energy'
        self.network = None

         
        
    def evaluate_thrust(self,state):
        
        """Evaluate the thrust produced by the energy network.
    
                Assumptions:
                Network has an "evaluate_thrust" method.
    
                Source:
                N/A
    
                Inputs:
                State data container
    
                Outputs:
                Results of the thrust evaluation method.
    
                Properties Used:
                N/A                
            """
                
            
        network = self.network
        results = network.evaluate_thrust(state) 
        
        return results


    def initialize(self):
        """Initializes the surrogate needed for the energy network calculations.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        self.geometry
        """  
        '''try:
            surrogate = self.network.turbofan
        except:
            try:
                surrogate = self.network.turbofan_hybrid
            except:
                surrogate = self.network.network 

       if surrogate.use_surrogate is True :
            efficiency_map = 'prop_eff.dat'

            # read the efficiency map
            with open(efficiency_map) as textFile:
                lines = [line.split() for line in textFile]
                lines = np.asarray(lines)
            lines = lines.astype(np.float)

            # split data into the grid and results
            xy       = lines[:,0:2]
            eta_data = lines[:,2:3] 

            # normalize the altitude
            max_alt = np.amax(xy[:,0])
            xy[:,0] = xy[:,0]/max_alt
    
            # set the surrogate data limits
            xlimits       = np.array([[np.amin(xy[:,0]),np.amax(xy[:,0])+0.1],[0, np.amax(xy[:,1])+0.1]])
            eta_surrogate = RMTB(print_global=False, num_ctrl_pts=60, xlimits=xlimits, nonlinear_maxiter=100, energy_weight=1e-12)
            #eta_surrogate.set_training_values(xy, np.transpose([eta_data]))
            eta_surrogate.set_training_values(xy, eta_data)
            eta_surrogate.train()

            self.network.network.eta_surrogate = Data()
            self.network.network.eta_surrogate.fit     = eta_surrogate
            self.network.network.eta_surrogate.max_alt = max_alt

            # plot surrogate data
            points_alt  = np.linspace(0.,np.amax(xy[:,0]),100)
            points_mach = np.linspace(0.,np.amax(xy[:,1]),100)

            num_alt  = len(points_alt)
            num_mach = len(points_mach)

            alt_mesh,mach_mesh = np.meshgrid(points_alt,points_mach)

            eta_sur = np.zeros(np.shape(alt_mesh))
                            
            x = np.zeros((num_alt,num_mach, 2))
            x[:, :, 0] = np.outer(np.linspace(0, 7800/max_alt, num_alt), np.ones(num_mach))
            x[:, :, 1] = np.outer(np.ones(num_alt), np.linspace(0.0, np.amax(xy[:,1]), num_mach))

            eta_sur = eta_surrogate.predict_values(x.reshape((num_alt * num_mach, 2)))[:, 0].reshape((num_alt, num_mach))
    
            plt.clf()
            #plt.plot( xy[:,0]*max_alt, xy[:,1], "o",)
            plt.contour(x[:, :, 0]*max_alt, x[:, :, 1], eta_sur, 20)
            plt.pcolormesh(x[:, :, 0]*max_alt, x[:, :, 1], eta_sur, cmap=plt.get_cmap("rainbow"))
            plt.xlabel("altitude, m")
            plt.ylabel("Mach number")
            plt.title("Propulsive efficiency surrogate")
            plt.colorbar()
            plt.savefig('surrogate_eta.png') '''
        
    finalize = initialize