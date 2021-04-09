## @ingroup Components-Energy-Networks
# Turboprop_LoFi.py
# 
# Created:  Feb 2021, S. Karpuk
# Modified: 
#           

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import numpy as np
from SUAVE.Components.Propulsors.Propulsor import Propulsor
from SUAVE.Core import Data, Units

# ----------------------------------------------------------------------
#  Network
# ----------------------------------------------------------------------
## @ingroup Components-Energy-Networks
class Turboprop_LoFi(Propulsor):
    """ A simple mock up of an internal combustion propeller engine. This network uses an imported 
        propeller efficiency map.
    
        Assumptions:
        None
        
        Source:
        None
    """      
    def __defaults__(self):
        """ This sets the default values for the network to function.
    
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
        self.engine            = None
        self.propeller         = None
        self.engine_length     = None
        self.number_of_engines = None
        self.thrust_angle      = 0.0
    
    # manage process with a driver function
    def evaluate_thrust(self,state):
        """ Calculate thrust given the current state of the vehicle
    
            Assumptions:
    
            Source:
            N/A
    
            Inputs:
            state [state()]
    
            Outputs:
            results.thrust_force_vector [newtons]
            results.vehicle_mass_rate   [kg/s]
            conditions.propulsion:
                power                [W]
    
            Properties Used:
            Defaulted values
        """           
        # unpack
        conditions  = state.conditions
        numerics    = state.numerics
        engine      = self.engine
        propeller   = self.propeller
        gearbox     = self.gearbox
        
        # Gearbox efficiency
        eta_g = gearbox.efficiency

        # Throttle the engine
        eta = conditions.propulsion.throttle[:,0,None]
        conditions.propulsion.combustion_engine_throttle = eta          # keep this 'throttle' on
        
        # Free-strem conditiona
        altitude = conditions.freestream.altitude     
        M        = conditions.freestream.mach_number 
        u0       = conditions.freestream.velocity

        # Propeller properties
        eta_surrogate = propeller.efficiency_map
        max_alt       = propeller.max_altitude

        # Run the engine
        engine.power(conditions)
        power_output     = engine.outputs.power
        max_power_output = engine.outputs.maximum_power  
        sigma            = engine.outputs.sigma
        mdot             = engine.outputs.fuel_flow_rate 
        
        # predict the propeller efficiency
        x    = np.concatenate((altitude/max_alt, M), axis=1)
        etap = eta_surrogate.predict_values(x) 

        F = np.divide(np.multiply(eta_g*etap,power_output),u0)
 
        conditions.propulsion.power                = power_output
        conditions.propulsion.maximum_power        = max_power_output
        conditions.propulsion.power_drop           = sigma
        conditions.propulsion.propeller_efficiency = etap
        
        # Create the outputs
        F    = self.number_of_engines * F * [np.cos(self.thrust_angle),0,-np.sin(self.thrust_angle)]      
        mdot *= self.number_of_engines

        results = Data()
        results.thrust_force_vector = F
        results.maximum_power       = max_power_output
        results.power               = power_output
        results.vehicle_mass_rate   = mdot
        
        return results
    
    
