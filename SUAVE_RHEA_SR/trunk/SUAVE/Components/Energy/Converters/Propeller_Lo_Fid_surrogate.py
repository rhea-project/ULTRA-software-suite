## @ingroup Components-Energy-Converters
# Propeller_Lo_Fid_surrogate.py
#
# Created:  Jun 2020, E. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import numpy as np
from SUAVE.Components.Energy.Energy_Component import Energy_Component
from warnings import warn
from SUAVE.Core import Data, Units
from smt.surrogate_models import RMTB

# ----------------------------------------------------------------------
#  Propeller Class
# ----------------------------------------------------------------------    
## @ingroup Components-Energy-Converters
class Propeller_Lo_Fid_surrogate(Energy_Component):
    """This is a low-fidelity propeller component using a surrogate dataset.
    
    Assumptions:
    None

    Source:
    None
    """    
    def __defaults__(self):
        """This sets the default values for the component to function.

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
        self.tip_radius            = 0.0
        self.propulsive_efficiency = 0.0

   
    def spin(self,conditions,propeller,motor,eta_surrogate,max_alt):
        """Analyzes a propeller given geometry and operating conditions.

        Assumptions:
        per source

        Source:
        Qprop theory document

        Inputs:
        self.inputs.omega            [radian/s]
        self.inputs.torque           [Nm]
        conditions.freestream.
          density                    [kg/m^3]
          dynamic_viscosity          [kg/(m-s)]
          velocity                   [m/s]
          speed_of_sound             [m/s]
          temperature                [K]

        Outputs:
        conditions.propulsion.etap   [-]  (propulsive efficiency)
        thrust                       [N]
        power                        [W]

        Properties Used:
        self.tip_radius              [m]
        self.propulsive_efficiency   [-]
        """    
           
        # Unpack
        max_power   = motor.sea_level_power
        altitude    = conditions.freestream.altitude[:,0,None]      
        M           = conditions.freestream.mach_number[:,0,None]  
        V           = conditions.freestream.velocity[:,0,None]
        throttle    = conditions.propulsion.throttle 

        x    = np.concatenate((altitude/max_alt, M), axis=1)
        etap = eta_surrogate.predict_values(x)        


        # Condition for negative throttle
        #for i in range(np.shape(throttle)[0]):
        #    if throttle[i,0] < 0:
        #        throttle[i,0] = 0
      
        if motor.tag != 'cryo-cooled motor':   
            # Calculate density ratio for the cooling power lapse calculation
            sigma = (1-0.000022558*altitude)**4.2561
            # calculate the power drop with the altitude for air-cooled engines
            power  = (max_power*1e3) * (0.4505*sigma+0.5519) * throttle                 # Accounts for the altitude effect and cooling (sqrt(density ratio))
        else:
            power  = max_power*1e3  * throttle 

        thrust = np.divide(np.multiply(etap,power),V)
        
        
        conditions.propulsion.etap = etap

        # store data
        results_conditions = Data     
        outputs   = results_conditions(

            #speed_of_sound            = conditions.freestream.speed_of_sound,
            density                   = conditions.freestream.density,                     
            power                     = power,
            propeller_efficiency      = etap,
            thrust                    = thrust
        ) 
      
        return  outputs, etap


    
