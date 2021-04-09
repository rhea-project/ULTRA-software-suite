## @ingroup Methods-Performance
# estimate_cruise_performance.py
#
# Created:  Sept, 2020 
# Modified: 
#


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from   SUAVE.Core import Data, Units

from scipy import optimize

from SUAVE.Analyses.Atmospheric         import US_Standard_1976 as atmosphere
from SUAVE.Methods.Power.Battery.Sizing import initialize_from_mass


import numpy as np

# ----------------------------------------------------------------------
#  Compute Cruise performance
# ----------------------------------------------------------------------

## @ingroup Methods-Performance
def flight_envelope(analyses,vehicle,altitude):
    """ Creates a flight envelope

    Assumptions:
    Efficiency of 85% is assumed (Constant speed propeller)

    Source:
    Gudmundsson, S., "General Aviation AIrcraft Design. Applied Methods and Procedures", 2013 

    Inputs:
    airport.
      atmosphere                           [SUAVE data type]
      altitude                             [m]
      delta_isa                            [K]
    vehicle.
      mass_properties.landing              [kg]
      reference_area                       [m^2]
      maximum_lift_coefficient (optional)  [Unitless]

    Outputs:
    landing_field_length                   [m]

    Properties Used:
    N/A
    """       
   
    # ==============================================
    # Unpack
    # ==============================================


    # return
    return 


def compute_stall_speed(analyses,vehicle,altitude): 
    """ Computes cruise stall speed

    Assumptions:

    Source:
    Gudmundsson, S., "General Aviation AIrcraft Design. Applied Methods and Procedures", 2013 

    Inputs:
    airport.
      atmosphere                           [SUAVE data type]
      altitude                             [m]
      delta_isa                            [K]
    vehicle.
      mass_properties.landing              [kg]
      reference_area                       [m^2]
      maximum_lift_coefficient (optional)  [Unitless]

    Outputs:
    landing_field_length                   [m]

    Properties Used:
    N/A
    """    

    return


def compute_maximum_cruise_speed(analyses,vehicle,altitude): 
    """ Computes maximum cruise speed

    Assumptions:

    Source:
    Gudmundsson, S., "General Aviation AIrcraft Design. Applied Methods and Procedures", 2013 

    Inputs:
    airport.
      atmosphere                           [SUAVE data type]
      altitude                             [m]
      delta_isa                            [K]
    vehicle.
      mass_properties.landing              [kg]
      reference_area                       [m^2]
      maximum_lift_coefficient (optional)  [Unitless]

    Outputs:
    landing_field_length                   [m]

    Properties Used:
    N/A
    """    

    # Unpack inputs
    W0 = vehicle.mass_properties.max_takeoff 
    AR = vehicle.wings['main_wing'].aspect_ratio 
    S  = vehicle.wings['main_wing'].areas.reference

    # constants
    g = 9.81;       etap = 0.85 

    # obtain atmospheric data
    atmo      = atmosphere()
    atmo_data = atmo.compute_values(altitude)
    p  = atmo_data.pressure[0,0]
    T  = atmo_data.temperature[0,0]
    ro = atmo_data.density[0,0]
    a  = atmo_data.speed_of_sound[0,0]
    mu = atmo_data.dynamic_viscosity[0,0]

    # Initialize conditions
    state = SUAVE.Analyses.Mission.Segments.Conditions.State()
    conditions = state.conditions
    conditions.update( SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics())
    
    state.conditions.freestream.altitude    = np.array([[altitude]])
    state.conditions.freestream.density     = np.array([[ro]])
    conditions.freestream.gravity           = np.array([[g]])  
    conditions.freestream.speed_of_sound    = np.array([[a]])
    conditions.freestream.temperature       = np.array([[T]])
    conditions.freestream.pressure          = np.array([[p]])
    conditions.freestream.dynamic_viscosity = np.array([[mu]])
    conditions.propulsion.throttle          = np.array([[1]])

    time_list     = [0, 0.5]
    time_int_list = [0.5]

    state.numerics.time = time_list
    state.numerics.time = Data()
    state.numerics.time.integrate = time_int_list

    aerodynamics = analyses.configs.base.aerodynamics
    aerodynamics.geometry = vehicle
    aerodynamics.initialize()

    # Obtain maximum power and propeller efficiency at the altitude 
    results = vehicle.propulsors.evaluate_thrust(state)    
    Power   = results.power[0,0]        


    # Convert units from SI to US
    W0    = 2.20462 * W0
    Power = 0.00134102 * Power
    S     = 10.7639 * S
    ro    = 0.00194032 * ro

    # Solve the cubic equation using the bisection method (in US units)
    count = 0;      eps = 10e-4;    error = 1
    V0    = 0.1;    V1  = 1000;   
    F0    = -1100 * etap * Power

    while (abs(error) > eps) and (count < 100):
        count += 1
        Vmid = 0.5*(V0 + V1)

        # Calculate CDmin
        Vmid_SI = Vmid*Units['ft/s']
        ro_SI   = ro / 0.00194032
        conditions.freestream.velocity          = np.array([[Vmid_SI]])
        conditions.freestream.mach_number       = np.array([[Vmid_SI/a]])
        conditions.freestream.reynolds_number   = np.array([[ro_SI*Vmid_SI/mu]])
        conditions.freestream.dynamic_pressure  = 0.5 * ro_SI * Vmid_SI**2

        aero_flight = aerodynamics.evaluate(state)
        CDmin       = aero_flight.drag.parasite.total[0,0] + aero_flight.drag.miscellaneous[0,0]
        e           = state.conditions.aerodynamics.drag_breakdown.induced.efficiency_factor[0,0]
        k           = 1/(np.pi*e*AR)

        # Calculate the function of the mid-point
        K1 = (550*etap*Power)**2 - (2*W0*Vmid)**2*CDmin*k
        if K1 < 0:
            K1 = 0
        Fmid = ro*S*CDmin*Vmid**3 - 550*etap*Power - np.sqrt(K1)

        if F0*Fmid < 0:
            V1 = Vmid
        else:
            V0 = Vmid
            F0 = Fmid

        error = (V1 - V0)/V0

    Vmax = 0.3048 * Vmid
    Mmax = Vmax / a

    # Write output
    fid = open('Max_cruise_speed.dat','w')   # Open output file
    fid.write('Output file with MAximum cruise speed at an altitude\n\n') #Start output printing
    fid.write( ' Cruise Altitude ...............: ' + str( '%8.2F'   %   altitude   ) + ' kg\n' )
    fid.write( ' Maximum cruise speed ..........: ' + str( '%8.2F'   %   Vmax    ) + ' m/s\n' )
    fid.write( ' Maximum Mach number ...........: ' + str( '%8.3F'   %   Mmax  ) + ' \n' )
    fid.close()

    return Vmax, Mmax