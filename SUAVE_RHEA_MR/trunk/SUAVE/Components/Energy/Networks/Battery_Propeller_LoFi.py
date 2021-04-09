## @ingroup Components-Energy-Networks
# Battery_Propeller_LoFi.py
# 
# Created:  Jun 2020, S.Karpuk
# Modified: 
#  
# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import scipy
import numpy as np
from SUAVE.Components.Propulsors.Propulsor import Propulsor
from smt.surrogate_models import RMTB
from SUAVE.Core import Data , Units

# ----------------------------------------------------------------------
#  Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Battery_Propeller_LoFi(Propulsor):
    """ This is a simple network with a battery powering a propeller through
        an electric motor
    
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
        self.propulsor        = None
        self.battery          = None 
        self.esc              = None
        self.inverter         = None
        self.tag              = 'Network'
        self.thrust_angle      = 0.0
        self.use_surrogate     = False
  


    # manage process with a driver function
    def evaluate_thrust(self,state):
        """ Calculate thrust given the current state of the vehicle
    
            Assumptions:
            Caps the throttle at 110% and linearly interpolates thrust off that
    
            Source:
            N/A
    
            Inputs:
            state [state()]
    
            Outputs:
            results.thrust_force_vector [newtons]
            results.vehicle_mass_rate   [kg/s]
            conditions.propulsion:
                rpm                  [radians/sec]
                current              [amps]
                battery_draw         [watts]
                battery_energy       [joules]
                voltage_open_circuit [volts]
                voltage_under_load   [volts]
                motor_torque         [N-M]
                propeller_torque     [N-M]
    
            Properties Used:
            Defaulted values
        """          
 
        # unpack
        conditions  = state.conditions
        numerics    = state.numerics
        motor       = self.motor
        propeller   = self.propeller
        thrust      = self.thrust
        esc         = self.esc
        inv         = self.inverter
        battery     = self.battery
        num_engines = self.number_of_engines
        I           = numerics.time.integrate
     

        # Set battery energy
        battery.current_energy = conditions.propulsion.battery_energy
        
        # Calculate thrust
        results, etap = propeller.spin(conditions,propeller,thrust,self.eta_surrogate.fit ,self.eta_surrogate.max_alt)

        # Power with the motor       
        propulsive_power    = np.reshape(results.power, (-1,1))
        motor_power         = propulsive_power/motor.motor_efficiency
        
        # Power with ESC
        ESC_power = motor_power/esc.efficiency
        
        # Power with Inverter
        Inv_power = ESC_power/inv.efficiency

        # Calculate battery draw and remaining battery energy
        #battery_draw   = scipy.integrate.cumtrapz(conditions.frames.inertial.time,ESC_power, initial = 0)     

        battery_draw = np.dot(I,Inv_power)

        # Pack the conditions for outputs 
        battery_energy = battery.current_energy 
        conditions.propulsion.battery_draw   = battery_draw
        conditions.propulsion.battery_energy = battery_energy - battery_draw

        # Create the outputs
        F    = results.thrust * [np.cos(self.thrust_angle),0,-np.sin(self.thrust_angle)]      
        mdot = np.zeros_like(F)

        F_mag = np.atleast_2d(np.linalg.norm(F, axis=1)/Units.lbs) # lb                       
        conditions.propulsion.power_loading         = (F_mag.T)/(battery_draw/Units.hp)                # lb/hp 
        
        results = Data()
        results.thrust_force_vector  = F
        results.vehicle_mass_rate    = mdot
        results.propeller_efficiency = etap 
        results.power                = propulsive_power
        
        return results
      
            
    __call__ = evaluate_thrust


