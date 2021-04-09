## @ingroup Components-Energy-Networks
#Turbofan.py
# 
# Created:  Oct 2014, A. Variyar, 
# Modified: Feb 2016, M. Vegh
#           Jul 2017, M. Clarke
#           Aug 2017, E. Botero
#           Oct 2017, E. Botero
#           Nov 2018, T. MacDonald

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import os

# package imports
import numpy as np

from SUAVE.Core import Data
from SUAVE.Components.Propulsors.Propulsor import Propulsor

# ----------------------------------------------------------------------
#  Turbofan Network
# ----------------------------------------------------------------------
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

## @ingroup Components-Energy-Networks
class Turbofan_Hybrid(Propulsor):
    """ This is a turbofan. 
    
        Assumptions:
        None
        
        Source:
        Most of the componentes come from this book:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/
    """      
    
    def __defaults__(self, model=None):
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
        
        #setting the default values
        self.tag = 'Turbofan'
        self.number_of_engines    = 1.0
        self.nacelle_diameter     = 1.0
        self.engine_length        = 1.0
        self.bypass_ratio         = 1.0
        self.SFC_adjustment       = 0.0 # Less than 1 is a reduction
        self.OpenVSP_flow_through = False
        
        #areas needed for drag; not in there yet
        self.areas             = Data()
        self.areas.wetted      = 0.0
        self.areas.maximum     = 0.0
        self.areas.exit        = 0.0
        self.areas.inflow      = 0.0
        self.sealevel_static_thrust = 0
        
        self.throttle_turbofan = []
        self.throttle_ductedfan = []        
        
        self.model = model
        
    _component_root_map = None
        
    # linking the different network components
    def evaluate_thrust(self, state):
        """ Calculate thrust given the current state of the vehicle
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            state [state()]
    
            Outputs:
            results.thrust_force_vector [newtons]
            results.vehicle_mass_rate   [kg/s]
            conditions.propulsion.acoustic_outputs:
                core:
                    exit_static_temperature      
                    exit_static_pressure       
                    exit_stagnation_temperature 
                    exit_stagnation_pressure
                    exit_velocity 
                fan:
                    exit_static_temperature      
                    exit_static_pressure       
                    exit_stagnation_temperature 
                    exit_stagnation_pressure
                    exit_velocity 
    
            Properties Used:
            Defaulted values
        """
        
         
        # unpack
        conditions = state.conditions
        
        hybridModel = self.model
        
        hybridModel.input_sim_vars = list(hybridModel.mod.getInputs().keys()) 
        dummy = hybridModel.input_sim_vars[0]
        
        # set input for hybridModel model
        hybridModel.profile.thrust = np.asarray(conditions.propulsion.throttle)
        hybridModel.profile.altitude = np.asarray(conditions.freestream.altitude)
        hybridModel.profile.MachNumber = np.asarray(conditions.freestream.mach_number)
        
        # get SUAVE time vector
        time_set = conditions.frames.inertial.time # SUAVE time vector
        start_time = time_set[0] # SUAVE start time
        final_time = time_set[-1] # SUAVE start time
        
        # if one point in time, generate tuple of input for hybridModel model
        tag = 0        
        if len(time_set) == 1:
            tag = 1
            hybridModel.profile.time = np.append(time_set, time_set+1) - start_time
            
            for i in hybridModel.input_sim_vars:
                hybridModel.profile[i] = np.append(hybridModel.profile[i], 
                                                   hybridModel.profile[i])
         
        # adjust time vector, so that start time is zero
        hybridModel.profile.time = time_set - start_time
        
        # set start and final time for hybridModel model                
        hybridModel.input_object = hybridModel.convertInput()
        
        start_time_temp = hybridModel.input_object[dummy][0][0]
        final_time_temp = hybridModel.input_object[dummy][-1][0]

        hybridModel.setInput()
        
        # set simulation options only, if start and stop time changes
        if isclose(start_time, hybridModel.start_time_set) & isclose(final_time, hybridModel.final_time_set):
            pass
        else:
            hybridModel.segment = hybridModel.segment + 1
            
            hybridModel.start_time_set = start_time
            hybridModel.final_time_set = final_time
            
            hybridModel.SOE_Jet_start = hybridModel.SOE_Jet_end
            hybridModel.SOE_Bat_start = hybridModel.SOE_Bat_end
            
            hybridModel.mod.setParameters(**{'Aircraft_Model.keroseneTankStateOfEnergyInit': hybridModel.SOE_Jet_start})
            hybridModel.mod.setParameters(**{'Aircraft_Model.Battery_NU.stateOfEnergyInit': hybridModel.SOE_Bat_start})
            
            hybridModel.start_time = start_time_temp
            hybridModel.final_time = final_time_temp
            step_size = (hybridModel.final_time - hybridModel.start_time) / 100
            hybridModel.mod.setSimulationOptions(startTime=hybridModel.start_time, 
                                                  stopTime=hybridModel.final_time, 
                                                  stepSize=step_size, 
                                                  tolerance=1e-06)

        # simulate:
        hybridModel.simulate()
        
        # get results
        res = hybridModel.read_mat(hybridModel.RES_NAME)
        
        time = res['time'] + start_time # reset time vector back to SUAVE start time
        thrustOperation = res['Aircraft_Model.TMS.thrustOperation']
        massFlowKerosene = res['Aircraft_Model.EMS_Kerosene.massFlowFuelSet'] 
        throttleTurbofan = res['Aircraft_Model.Gas_Turbine.Operation.throttleOperation'] 
        throttleDuctedfan = res['Aircraft_Model.Ducted_Fan_Elec.Operation.throttleOperation'] 
        thrust = thrustOperation
                
        hybridModel.SOE_Jet_end = res['Aircraft_Model.Kerosene_Tank_PU.CMS.stateOfEnergy'][-1]
        hybridModel.SOE_Bat_end = res['Aircraft_Model.Battery_NU.Battery_PU.CMS.stateOfEnergy'][-1]
        
        electric_motor_power = res['Aircraft_Model.Electric_Motor_NU.Electric_Motor_PU.powerSet'][:-1]
        gas_turbine_power = res['Aircraft_Model.Gas_Turbine.powerSet'][:-1]
        battery_power = res['Aircraft_Model.Battery_NU.Battery_PU.powerSet'][:-1]
        kerosene_power = res['Aircraft_Model.Kerosene_Tank_PU.powerSet'][:-1]
        time_steps = np.diff(res['time'])
        
        hybridModel.data.thrust[hybridModel.segment] = thrustOperation        
        hybridModel.data.max_throttle_turbofan[hybridModel.segment] = max(throttleTurbofan)
        hybridModel.data.max_throttle_ductedfan[hybridModel.segment] = max(throttleDuctedfan)
        hybridModel.data.max_power_motor[hybridModel.segment] = max(electric_motor_power)
        hybridModel.data.max_power_battery[hybridModel.segment] = max(battery_power)
        hybridModel.data.max_power_gas_turbine[hybridModel.segment] = max(gas_turbine_power)
        
        hybridModel.data.min_SOC_battery[hybridModel.segment] = hybridModel.SOE_Bat_end
        hybridModel.data.min_SOC_kerosene_tank[hybridModel.segment] = hybridModel.SOE_Jet_end
        hybridModel.data.battery_energy[hybridModel.segment] = np.dot(time_steps, battery_power) / 3600
        hybridModel.data.kerosene_energy[hybridModel.segment] = np.dot(time_steps, kerosene_power) / 3600
                
        if tag == 1:
            thrust = np.array([thrust[0]])
            massFlowKerosene = np.array([massFlowKerosene[0]])
            throttle = np.array([throttleTurbofan[0]])
        else:
            # interpolate results on SUAVE time
            thrust = np.interp(time_set, time, thrust)
            massFlowKerosene = np.interp(time_set, time, massFlowKerosene)
            throttle = np.interp(time_set, time, throttleTurbofan)
        
        # create numpy ndarray for results
        thrust_nd = np.ndarray(shape=(len(thrust), 1))
        massFlowKerosene_nd = np.ndarray(shape=(len(massFlowKerosene), 1))
        throttle_nd = np.ndarray(shape=(len(throttle), 1))
        
        for i in range(len(thrust_nd)):
            thrust_nd[i] = thrust[i]
            massFlowKerosene_nd[i] = massFlowKerosene[i]
            throttle_nd[i] = throttle[i]
    
        # getting the network outputs from the thrust outputs        
        F            = thrust_nd*[1,0,0]
        mdot         = massFlowKerosene_nd
        F_vec        = conditions.ones_row(3) * 0.0
        F_vec[:,0]   = F[:,0]
        F            = F_vec

        results = Data()
        results.thrust_force_vector = F
        results.vehicle_mass_rate   = mdot
        
        conditions.propulsion.throttle = throttle_nd
		
        return results
    
    def size(self,state):  
        """ Size the turbofan
    
            Assumptions:
            None
    
            Source:
            N/A
    
            Inputs:
            State [state()]
    
            Outputs:
            None
    
            Properties Used:
            N/A
        """             
        
        #Unpack components
        conditions = state.conditions
        thrust     = self.thrust
        thrust.size(conditions)
        
    def engine_out(self,state):
        """ Lose an engine
    
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
        
        temp_throttle = np.zeros(len(state.conditions.propulsion.throttle))
        
        for i in range(0,len(state.conditions.propulsion.throttle)):
            temp_throttle[i] = state.conditions.propulsion.throttle[i]
            state.conditions.propulsion.throttle[i] = 1.0
        
        results = self.evaluate_thrust(state)
        
        for i in range(0,len(state.conditions.propulsion.throttle)):
            state.conditions.propulsion.throttle[i] = temp_throttle[i]
        
        results.thrust_force_vector = results.thrust_force_vector/self.number_of_engines*(self.number_of_engines-1)
        results.vehicle_mass_rate   = results.vehicle_mass_rate/self.number_of_engines*(self.number_of_engines-1)

        return results
        
    __call__ = evaluate_thrust