# Turboprop_LoFi.py
#
# Created:  Feb 2021: S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
import SUAVE
from SUAVE.Core import Data, Units

# package imports
import numpy as np
from SUAVE.Components.Energy.Energy_Component import Energy_Component

# ----------------------------------------------------------------------
#  Internal Combustion Engine Class
# ----------------------------------------------------------------------

class Turboprop_LoFi(Energy_Component):

    def __defaults__(self):

        self.sea_level_power    = 0.0
        self.flat_rate_altitude = 0.0
        self.speed              = 0.0
        self.power_specific_fuel_consumption = 0.5 # lb/hr/hp :: Ref: Table 3.4, Daniel Raymer, Aircraft Design: Conceptual Approach

    def power(self,conditions):
        """ The turboprop engine output power and specific power consumption
        Inputs:
            Engine:
                sea-level power
                flat rate altitude
                throttle setting
            Freestream conditions:
                altitude
                delta_isa
        Outputs:
            Brake power (or Shaft power)
            Power (brake) specific fuel consumption
            Fuel flow
            Torque
        """

        # Unpack
        altitude  = conditions.freestream.altitude
        delta_isa = conditions.freestream.delta_ISA
        throttle  = conditions.propulsion.combustion_engine_throttle
        PSLS      = self.sea_level_power
        h_flat    = self.flat_rate_altitude
        speed     = self.speed
        BSFC      = self.power_specific_fuel_consumption
        n         = 0.75            # a power reduction factor from the Ruijgrok method
        #n         = 0.7             # a power reduction factor from the Anderson's method

        altitude_virtual = altitude - h_flat # shift in power lapse due to flat rate
        atmo = SUAVE.Analyses.Atmospheric.US_Standard_1976()
        atmo_values = atmo.compute_values(altitude_virtual,delta_isa)
        #atmo_values = atmo.compute_values(altitude_virtual)
        p   = atmo_values.pressure
        T   = atmo_values.temperature
        rho = atmo_values.density
        a   = atmo_values.speed_of_sound
        mu  = atmo_values.dynamic_viscosity

        # computing the sea-level ISA atmosphere conditions
        atmo_values = atmo.compute_values(0,0)
        p0   = atmo_values.pressure[0,0]
        T0   = atmo_values.temperature[0,0]
        rho0 = atmo_values.density[0,0]
        a0   = atmo_values.speed_of_sound[0,0]
        mu0  = atmo_values.dynamic_viscosity[0,0]

        # calculating the density ratio:
        sigma = rho / rho0

        # calculating available power based on Gagg and Ferrar model (ref: S. Gudmundsson, 2014 - eq. 7-16)
        Pavailable = PSLS * 1e3 * sigma**n      
        Pavailable[h_flat > altitude] = PSLS * 1e3

        # applying throttle setting
        output_power = Pavailable * throttle

        output_power[output_power<0.] = 0.

        SFC = BSFC * Units['lb/hp/hr']

        #fuel flow rate
        a = np.zeros_like(altitude)
        fuel_flow_rate   = np.fmax(output_power*SFC,a)

        # store to outputs
        self.outputs.power                           = output_power
        self.outputs.maximum_power                   = Pavailable
        self.outputs.sigma                           = sigma**n
        self.outputs.power_specific_fuel_consumption = BSFC
        self.outputs.fuel_flow_rate                  = fuel_flow_rate

        return self.outputs
