## @ingroup Components-Energy-Networks
#Turbofan.py
# 
# Created:  Oct, S. Karpuk, 
# Modified: 
#


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE

# package imports
import numpy as np

from SUAVE.Core import Data
from SUAVE.Components.Propulsors.Propulsor import Propulsor

# ----------------------------------------------------------------------
#  Turbofan Network
# ----------------------------------------------------------------------

## @ingroup Components-Energy-Networks
class Turbofan_DuctedFan_LoFi(Propulsor):
    """ Low-fidelity parallel Turbofan - Ducted Fan propulsion system model. 
    
        Assumptions:
        None
        
        Source:
        Most of the turbofan componentes come from this book:
        https://web.stanford.edu/~cantwell/AA283_Course_Material/AA283_Course_Notes/

        The electric ducted fan:



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
        
        #setting the default values
        self.tag = 'Hybrid_Turbofan_DuctedFan'
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
    _component_root_map = None
        
    # linking the different network components
    def evaluate_thrust(self,state):
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

        #Unpack
        #------------------------------------------------
        conditions = state.conditions
        
        ram                        = self.ram
        # Unpack turbofan components
        inlet_nozzle               = self.inlet_nozzle
        low_pressure_compressor    = self.low_pressure_compressor
        high_pressure_compressor   = self.high_pressure_compressor
        fan                        = self.fan
        combustor                  = self.combustor
        high_pressure_turbine      = self.high_pressure_turbine
        low_pressure_turbine       = self.low_pressure_turbine
        core_nozzle                = self.core_nozzle
        fan_nozzle                 = self.fan_nozzle
        thrust                     = self.thrust
        bypass_ratio               = self.bypass_ratio

        num_engines = self.number_of_engines

        inlet_nozzle_DF              = self.inlet_nozzle_DF
        fan_DF                       = self.fan_DF
        fan_nozzle_DF                = self.fan_nozzle_DF
        thrust_DF                    = self.thrust_DF
        
        Hp     = thrust.hybridization
        phi    = thrust.control_parameter
        Tmax   = thrust.total_design

        # Split max thrust based on the degree of hybridization
        TmaxGT = Tmax * (1-Hp)
        TmaxDF = Tmax * Hp 


        #Creating the network by manually linking the different components
        
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
        ram(conditions) 

        
        # Calculate total thrust assuming all engines as turbofans
        #---------------------------------------------------------------------------------------------------
        #link inlet nozzle to ram 
        inlet_nozzle.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
        inlet_nozzle.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure
        
        #Flow through the inlet nozzle
        inlet_nozzle(conditions)

        #--link low pressure compressor to the inlet nozzle
        low_pressure_compressor.inputs.stagnation_temperature  = inlet_nozzle.outputs.stagnation_temperature
        low_pressure_compressor.inputs.stagnation_pressure     = inlet_nozzle.outputs.stagnation_pressure
        
        #Flow through the low pressure compressor
        low_pressure_compressor(conditions)

        #link the high pressure compressor to the low pressure compressor
        high_pressure_compressor.inputs.stagnation_temperature = low_pressure_compressor.outputs.stagnation_temperature
        high_pressure_compressor.inputs.stagnation_pressure    = low_pressure_compressor.outputs.stagnation_pressure
        
        #Flow through the high pressure compressor
        high_pressure_compressor(conditions)
        
        #Link the fan to the inlet nozzle
        fan.inputs.stagnation_temperature                      = inlet_nozzle.outputs.stagnation_temperature
        fan.inputs.stagnation_pressure                         = inlet_nozzle.outputs.stagnation_pressure
        
        #flow through the fan
        fan(conditions)
        
        #link the combustor to the high pressure compressor
        combustor.inputs.stagnation_temperature                = high_pressure_compressor.outputs.stagnation_temperature
        combustor.inputs.stagnation_pressure                   = high_pressure_compressor.outputs.stagnation_pressure
        
        #flow through the high pressor comprresor
        combustor(conditions)

        # link the shaft power output to the low pressure compressor
        try:
            shaft_power = self.Shaft_Power_Off_Take       
            shaft_power.inputs.mdhc                            = thrust.compressor_nondimensional_massflow
            shaft_power.inputs.Tref                            = thrust.reference_temperature
            shaft_power.inputs.Pref                            = thrust.reference_pressure
            shaft_power.inputs.total_temperature_reference     = low_pressure_compressor.outputs.stagnation_temperature
            shaft_power.inputs.total_pressure_reference        = low_pressure_compressor.outputs.stagnation_pressure
    
            shaft_power(conditions)
        except:
            pass

        #link the high pressure turbine to the combustor
        high_pressure_turbine.inputs.stagnation_temperature    = combustor.outputs.stagnation_temperature
        high_pressure_turbine.inputs.stagnation_pressure       = combustor.outputs.stagnation_pressure
        high_pressure_turbine.inputs.fuel_to_air_ratio         = combustor.outputs.fuel_to_air_ratio
        
        #link the high pressuer turbine to the high pressure compressor
        high_pressure_turbine.inputs.compressor                = high_pressure_compressor.outputs
        
        #link the high pressure turbine to the fan
        high_pressure_turbine.inputs.fan                       = fan.outputs
        high_pressure_turbine.inputs.bypass_ratio              = 0.0 #set to zero to ensure that fan not linked here
        
        #flow through the high pressure turbine
        high_pressure_turbine(conditions)
                
        #link the low pressure turbine to the high pressure turbine
        low_pressure_turbine.inputs.stagnation_temperature     = high_pressure_turbine.outputs.stagnation_temperature
        low_pressure_turbine.inputs.stagnation_pressure        = high_pressure_turbine.outputs.stagnation_pressure
        
        #link the low pressure turbine to the low_pressure_compresor
        low_pressure_turbine.inputs.compressor                 = low_pressure_compressor.outputs
        
        #link the low pressure turbine to the combustor
        low_pressure_turbine.inputs.fuel_to_air_ratio          = combustor.outputs.fuel_to_air_ratio
        
        #link the low pressure turbine to the fan
        low_pressure_turbine.inputs.fan                        = fan.outputs
        
        # link the low pressure turbine to the shaft power, if needed
        try:
            low_pressure_turbine.inputs.shaft_power_off_take   = shaft_power.outputs
        except:
            pass
        
        #get the bypass ratio from the thrust component
        low_pressure_turbine.inputs.bypass_ratio               = bypass_ratio
        
        #flow through the low pressure turbine
        low_pressure_turbine(conditions)
        
        #link the core nozzle to the low pressure turbine
        core_nozzle.inputs.stagnation_temperature              = low_pressure_turbine.outputs.stagnation_temperature
        core_nozzle.inputs.stagnation_pressure                 = low_pressure_turbine.outputs.stagnation_pressure
        
        #flow through the core nozzle
        core_nozzle(conditions)

        #link the dan nozzle to the fan
        fan_nozzle.inputs.stagnation_temperature               = fan.outputs.stagnation_temperature
        fan_nozzle.inputs.stagnation_pressure                  = fan.outputs.stagnation_pressure
        
        # flow through the fan nozzle
        fan_nozzle(conditions)
        
        # compute the thrust using the thrust component
        #link the thrust component to the fan nozzle
        thrust.inputs.fan_exit_velocity                        = fan_nozzle.outputs.velocity
        thrust.inputs.fan_area_ratio                           = fan_nozzle.outputs.area_ratio
        thrust.inputs.fan_nozzle                               = fan_nozzle.outputs
        
        #link the thrust component to the core nozzle
        thrust.inputs.core_exit_velocity                       = core_nozzle.outputs.velocity
        thrust.inputs.core_area_ratio                          = core_nozzle.outputs.area_ratio
        thrust.inputs.core_nozzle                              = core_nozzle.outputs
        
        #link the thrust component to the combustor
        thrust.inputs.fuel_to_air_ratio                        = combustor.outputs.fuel_to_air_ratio
        
        #link the thrust component to the low pressure compressor 
        thrust.inputs.total_temperature_reference              = low_pressure_compressor.outputs.stagnation_temperature
        thrust.inputs.total_pressure_reference                 = low_pressure_compressor.outputs.stagnation_pressure
        thrust.inputs.number_of_engines                        = num_engines
        thrust.inputs.bypass_ratio                             = bypass_ratio
        thrust.inputs.flow_through_core                        = 1./(1.+bypass_ratio) #scaled constant to turn on core thrust computation
        thrust.inputs.flow_through_fan                         = bypass_ratio/(1.+bypass_ratio) #scaled constant to turn on fan thrust computation        

        #compute the thrust
        thrust(conditions)

        # total thrust, power, and flow rate
        Ftot          = thrust.outputs.thrust*[1,0,0]
        mdot          = thrust.outputs.fuel_flow_rate
        F_vectot      = conditions.ones_row(3) * 0.0
        F_vectot[:,0] = Ftot[:,0]
        Ftot          = F_vectot
        #---------------------------------------------------------------------------------------------------

        # Calculate thrust distribution
        #---------------------------------------------------------------------------------------------------
        # Split the thrust based on the control parameter
        T   = thrust.outputs.thrust[0]
        TDF = np.max(Ftot-TmaxGT,0) + phi*(np.min(TmaxDF-Ftot,TmaxDF)-np.max(Ftot-TmaxGT,0))
        TGT = T - TDF

        # Update fuel flow based on the gas turbine thrust
        thrust.outputs.fuel_flow_rate = thrust.outputs.fuel_flow_rate * TGT/Ftot


        # Calculate total thrust assuming all engines as ducted fans (used only to calculate the ducted fan SFC)
        #---------------------------------------------------------------------------------------------------
#       Creating the network by manually linking the different components
        #set the working fluid to determine the fluid properties
        ram.inputs.working_fluid                               = self.working_fluid
        
        #Flow through the ram , this computes the necessary flow quantities and stores it into conditions
        ram(conditions)
        
        #link inlet nozzle to ram 
        inlet_nozzle_DF.inputs.stagnation_temperature             = ram.outputs.stagnation_temperature 
        inlet_nozzle_DF.inputs.stagnation_pressure                = ram.outputs.stagnation_pressure 
        
        #Flow through the inlet nozzle
        inlet_nozzle_DF(conditions)
        
        #Link the fan to the inlet nozzle
        fan_DF.inputs.stagnation_temperature                      = inlet_nozzle_DF.outputs.stagnation_temperature
        fan_DF.inputs.stagnation_pressure                         = inlet_nozzle_DF.outputs.stagnation_pressure
        
        #flow through the fan
        fan_DF(conditions)

        #link the fan nozzle to the fan
        fan_nozzle_DF.inputs.stagnation_temperature               = fan.outputs.stagnation_temperature
        fan_nozzle_DF.inputs.stagnation_pressure                  = fan.outputs.stagnation_pressure
        thrust_DF.inputs.fuel_to_air_ratio                        = 0.
        
         # flow through the fan nozzle
        fan_nozzle_DF(conditions)
        
        # compute the thrust using the thrust component
        #link the thrust component to the fan nozzle
        thrust_DF.inputs.fan_exit_velocity                        = fan_nozzle_DF.outputs.velocity
        thrust_DF.inputs.fan_area_ratio                           = fan_nozzle_DF.outputs.area_ratio
        thrust_DF.inputs.fan_nozzle                               = fan_nozzle_DF.outputs
        thrust_DF.inputs.bypass_ratio                             = 0.0
        
        #compute the thrust
        thrust_DF(conditions)

        Ftot_DF          = thrust.outputs.thrust*[1,0,0]
        F_vectot_DF      = conditions.ones_row(3) * 0.0
        F_vectot_DF[:,0] = Ftot_DF[:,0]
        Ftot_DF          = F_vectot_DF

        mdot_DF          = thrust_DF.outputs.fuel_flow_rate * TDF/Ftot_DF
        #---------------------------------------------------------------------------------------------------

        results = Data()
        results.thrust_force_vector = Ftot
        results.vehicle_mass_rate   = mdot + mdot_DF
        
        # store data
        results_conditions = Data
        conditions.propulsion.acoustic_outputs.core = results_conditions(
        exit_static_temperature             = core_nozzle.outputs.static_temperature,
        exit_static_pressure                = core_nozzle.outputs.static_pressure,
        exit_stagnation_temperature         = core_nozzle.outputs.stagnation_temperature,
        exit_stagnation_pressure            = core_nozzle.outputs.static_pressure,
        exit_velocity                       = core_nozzle.outputs.velocity
        )
        
        conditions.propulsion.acoustic_outputs.fan = results_conditions(
        exit_static_temperature             = fan_nozzle.outputs.static_temperature,
        exit_static_pressure                = fan_nozzle.outputs.static_pressure,
        exit_stagnation_temperature         = fan_nozzle.outputs.stagnation_temperature,
        exit_stagnation_pressure            = fan_nozzle.outputs.static_pressure,
        exit_velocity                       = fan_nozzle.outputs.velocity
        )
        
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