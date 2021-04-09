## @ingroup Methods-Performance
# estimate_climb_performance.py
#
# Created:  Aug, 2020 
# Modified: 
#


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from   SUAVE.Core import Data, Units

from SUAVE.Analyses.Atmospheric         import US_Standard_1976 as atmosphere
from SUAVE.Methods.Power.Battery.Sizing import initialize_from_mass


import numpy as np

# ----------------------------------------------------------------------
#  Compute field length required for landing
# ----------------------------------------------------------------------

## @ingroup Methods-Performance
def estimate_climb_performance(analyses,vehicle,altitude):
    """ Computes climb performance curves for a given airplane

    Assumptions:
    Is used only with AVL for now

    Source:
    Torenbeek, E., "Advanced Aircraft Design", 2013 (equation 9.25)

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

    weight  = vehicle.mass_properties.max_takeoff
    Sw      = vehicle.reference_area
    
    # constants
    g = 9.81;   CL_lim = 1.5

    # define the Mach number and the ROC sweeps
    mach = np.linspace(0.03, 0.85, num=100)
    Vx   = np.zeros(len(mach))
    Vy   = np.zeros(len(mach))

    # obtain atmospheric data
    atmo      = atmosphere()
    atmo_data = atmo.compute_values(altitude)
    p  = atmo_data.pressure[0,0]
    T  = atmo_data.temperature[0,0]
    ro = atmo_data.density[0,0]
    a  = atmo_data.speed_of_sound[0,0]
    mu = atmo_data.dynamic_viscosity[0,0]

    atmo_SL      = atmosphere()
    atmo_data_SL = atmo_SL.compute_values(0)
    p0  = atmo_data_SL.pressure[0,0]
    ro0 = atmo_data_SL.density[0,0]

    # compute airspeeds
    Vx = mach * atmo_data.speed_of_sound[0,0]

    qc  = p * (np.power(1+0.2*mach,3.5) - 1)
    Vx_eas = Vx * (ro/ro0)**0.5
    Vx_cas = np.multiply(Vx_eas,(p0/p)**0.5*np.power(np.divide((np.power(qc/p0+1,0.286)-1),(np.power(qc/p+1,0.286)-1)),0.5)) 

    q   = 0.5*ro*np.power(Vx,2)

    CL = weight*g/(q*Sw)

    #state = Data()
    #state = SUAVE.Analyses.Mission.Segments.Conditions.State()
    #state.conditions = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    #state.conditions.freestream = Data()

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

    for i in range(len(mach)):
        conditions.freestream.velocity          = np.array([[Vx[i]]])
        conditions.freestream.mach_number       = np.array([[Vx[i]/a]])
        conditions.freestream.reynolds_number   = np.array([[ro*Vx[i]/mu]])
        conditions.freestream.dynamic_pressure  = 0.5 * ro * Vx[i]**2
    
        # estimate the lift-curve slope
        AoA0 = 0.0 * Units.degrees;    AoA4 = 4.0 * Units.degrees; 
        state.conditions.aerodynamics.angle_of_attack = np.array([[AoA0]])   
        aero_aoa0 = aerodynamics.evaluate(state)
        state.conditions.aerodynamics.angle_of_attack = np.array([[AoA4]])  
        aero_aoa4 = aerodynamics.evaluate(state)

        CLa = (aero_aoa4.lift.total - aero_aoa0.lift.total) / (AoA4 - AoA0)

        # Calculate the flight AoA
        AoA = np.array([[CL[i]/CLa + AoA0]])

        state.conditions.aerodynamics.angle_of_attack = AoA
        aero_flight = aerodynamics.evaluate(state)
        Drag        = aero_flight.drag.total[0,0] * q[i] * Sw

        # Evaluate maximum thrust at a given flight condition
        results = vehicle.propulsors.evaluate_thrust(state) # total thrust  
        thrust  = np.linalg.norm(results.thrust_force_vector[0])

        Vy[i] = Vx[i] * (thrust - Drag) / (weight*g)

    Vx_data = np.vstack((Vx, Vx_eas, Vx_cas)) 

    # limit the data range based on the lift coefficient
    i = 0
    while CL[i] > CL_lim:
      i += 1

    Vy      = Vy[i::]
    Vx_data = Vx_data[:,i::]

    # find max ROC for a given altitude
    Vy_max   = np.max(Vy)  
    max_ind  = Vy.argmax()
    Vx_max   = Vx_data[:,max_ind]
    V_cl_max = np.hstack((Vx_max, Vy_max))

    # return
    return Vy, Vx_data, V_cl_max

  
