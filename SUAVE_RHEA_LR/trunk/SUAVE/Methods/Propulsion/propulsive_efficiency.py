## @ingroup Methods-Propulsion
# propulsive_efficiency.py
# 
# Created:  Mar 2020, S. Karpuk
# Modified: 
#          

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
import pylab as plt
from SUAVE.Core import Units
from smt.surrogate_models import RMTB
from SUAVE.Core import Data

# ----------------------------------------------------------------------
#  Propeller Design
# ----------------------------------------------------------------------
    
def propulsive_efficiency(efficiency_map):
    """ The script creates an SMT-map of propulsive efficiency based on the input data.
          
          Inputs:
          Either design power or thrust
          prop_attributes.
            efficiency_data                  [Unitless array]
                    format:
                        [altitude,m      Mach     efficiency]
                        |                                   |
                        |                                   |
                        [                                   ]
           
          Outputs:
            eta_surrogate                   [Unitless surface fit]
            max_alt                          maximum altitude, m
            
          Assumptions/ Source:

    """    

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
    plt.plot( xy[:,0]*max_alt, xy[:,1], "o",)
    plt.contour(x[:, :, 0]*max_alt, x[:, :, 1], eta_sur, 20)
    plt.pcolormesh(x[:, :, 0]*max_alt, x[:, :, 1], eta_sur, cmap=plt.get_cmap("rainbow"))
    plt.xlabel("altitude, m")
    plt.ylabel("Mach number")
    plt.title("eta surrogate")
    plt.colorbar()

    plt.savefig('surrogate_eta.png')                             

    return eta_surrogate,max_alt
