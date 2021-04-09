## @ingroup Methods-Performance
# estimate_take_off_field_length.py
#
# Created:  Jun 2014, T. Orra, C. Ilario, Celso, 
# Modified: Apr 2015, M. Vegh 
#           Jan 2016, E. Botero
#           Mar 2020, M. Clarke
#           May 2020, E. Botero
#           Jul 2020, E. Botero 


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUave Imports
import SUAVE
from SUAVE.Core import Data, Units

from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Helper_Functions import windmilling_drag
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Helper_Functions import estimate_2ndseg_lift_drag_ratio
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Helper_Functions import asymmetry_drag
from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Lift import compute_max_lift_coeff

from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Drag import parasite_drag_wing
from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Drag        import compute_induced_drag_highlift, compute_parasite_drag_highlift, \
                                                                 compute_interference_drag_highlift, \
                                                                 compute_landing_gear_drag    

# package imports
import numpy as np

# ----------------------------------------------------------------------
#  Compute field length required for takeoff
# ----------------------------------------------------------------------

## @ingroup Methods-Performance
def estimate_take_off_field_length(vehicle,analyses,airport,compute_2nd_seg_climb = 0):
    """ Computes the takeoff field length for a given vehicle configuration in a given airport.
    Also optionally computes the second segment climb gradient.

    Assumptions:
    For second segment climb gradient:
        1.One engine inoperative
    2. Oswald efficiency does not change with flap settings


    Source:
    http://adg.stanford.edu/aa241/AircraftDesign.html

    Inputs:
    analyses.base.atmosphere               [SUAVE data type]
    airport.
      altitude                             [m]
      delta_isa                            [K]
    vehicle.
      mass_properties.takeoff              [kg]
      reference_area                       [m^2]
      V2_VS_ratio (optional)               [Unitless]
      maximum_lift_coefficient (optional)  [Unitless]
      propulsors.*.number_of_engines       [Unitless]

    Outputs:
    takeoff_field_length                   [m]

    Properties Used:
    N/A
    """        

    # Constants
    surf_mu = 0.05              # road skin friction coefficient
    h_obst  = 50 * Units.ft 

    # ==============================================
    # Unpack
    # ==============================================
    atmo            = analyses.configs.takeoff.atmosphere
    altitude        = airport.altitude * Units.ft
    delta_isa       = airport.delta_isa
    weight          = vehicle.mass_properties.takeoff
    reference_area  = vehicle.reference_area
    b               = vehicle.wings['main_wing'].spans.projected
    dihedral        = vehicle.wings['main_wing'].dihedral
    LG_origin       = vehicle.landing_gear.main_origin 
    wing_origin     = vehicle.wings['main_wing'].origin 

    try:
        V2_VS_ratio = vehicle.V2_VS_ratio
    except:
        V2_VS_ratio = 1.20

    dtime = 1.5                 # time step for take-off calculation

    # ==============================================
    # Computing atmospheric conditions
    # ==============================================
    atmo_values       = atmo.compute_values(altitude,delta_isa)
    conditions        = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    
    p   = atmo_values.pressure
    T   = atmo_values.temperature
    rho = atmo_values.density
    a   = atmo_values.speed_of_sound
    mu  = atmo_values.dynamic_viscosity
    sea_level_gravity = atmo.planet.sea_level_gravity
    
    # ==============================================
    # Determining vehicle maximum lift coefficient
    # ==============================================
    # Condition to CLmax calculation: 90KTAS @ airport
    state = Data()
    state.conditions = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    state.conditions.freestream = Data()
    state.conditions.freestream.density           = rho
    state.conditions.freestream.velocity          = 90. * Units.knots
    state.conditions.freestream.dynamic_viscosity = mu
    
    state.conditions.aerodynamics.angle_of_attack = 0.0 
    
    settings = analyses.configs.takeoff.aerodynamics.settings

    maximum_lift_coefficient, zero_aoa_lift_coef, alphad0, c_pr_c, flap_area, slat_area, Cla = compute_max_lift_coeff(state,settings,vehicle)

    # ==============================================
    # Estimate the ground effect factor
    # ==============================================
    h  = np.average([wing_origin[2]-LG_origin[2], wing_origin[2]-LG_origin[2]+b*np.tan(dihedral)])
    hb = h/b

    phi_data    = np.zeros(3)
    phi_data[0] = 1 - (1-1.32*hb)/(1.05+7.4*hb)              # Weisselberger estimation
    phi_data[1] = (16*hb)**2/(1+(16*hb)**2)                  # McCormick estimation
    phi_data[2] = 1 - 2/np.pi**2 * np.log(1+np.pi/(8*hb))    # Asselin estimation
    phi = np.average(phi_data)

    # ==============================================
    # Computing speeds (Vs, V_LOF, V2)]
    # ==============================================
    stall_speed = (2 * weight * sea_level_gravity / (rho * reference_area * maximum_lift_coefficient)) ** 0.5
    V_LOF       = 1.1 * stall_speed
    V2_speed    = 1.2 * stall_speed

    # ==============================================
    # Determining vehicle number of engines
    # ==============================================
    engine_number = 0.
    for propulsor in vehicle.propulsors : # may have than one propulsor
        engine_number += propulsor.number_of_engines
    if engine_number == 0:
        raise ValueError("No engine found in the vehicle")

    # ==============================================
    # Initial conditions (stationary aircraft)
    # ============================================== 
    ground_distance = 0
    j = 0
    airspeed = 0.1 
    state = SUAVE.Analyses.Mission.Segments.Conditions.State()
    conditions = state.conditions
    conditions.update( SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics())

    '''conditions.freestream.dynamic_pressure  = 0.5 * rho * airspeed**2
    conditions.freestream.gravity           = np.array([[sea_level_gravity]])
    conditions.freestream.velocity          = np.array([[airspeed]])
    conditions.freestream.mach_number       = np.array(airspeed/a)
    conditions.freestream.speed_of_sound    = np.array([[a]])
    conditions.freestream.temperature       = np.array([[T]])
    conditions.freestream.pressure          = np.array([[p]])
    conditions.freestream.dynamic_viscosity = np.array([[mu]])
    conditions.freestream.reynolds_number   = np.array([[rho*airspeed/mu]])
    conditions.propulsion.throttle          = np.array([[1]])'''

    conditions.freestream.dynamic_pressure  = 0.5 * rho * airspeed**2
    conditions.freestream.gravity           = np.array(sea_level_gravity)
    conditions.freestream.velocity          = np.array(airspeed)
    conditions.freestream.mach_number       = np.array(airspeed/a)
    conditions.freestream.speed_of_sound    = np.array(a)
    conditions.freestream.temperature       = np.array(T)
    conditions.freestream.pressure          = np.array(p)
    conditions.freestream.dynamic_viscosity = np.array(mu)
    conditions.freestream.reynolds_number   = np.array(rho*airspeed/mu)
    conditions.propulsion.throttle          = np.array(1)

    time = 0
    battery_draw = 0.0
    while airspeed <= V_LOF:

        time_list     = [0, dtime]
        time_int_list = [dtime]

        state.numerics.time = time_list
        state.numerics.time = Data()
        state.numerics.time.integrate = time_int_list

        airspeed_old = airspeed
        results      = vehicle.propulsors.evaluate_thrust(state) # total thrust  
        thrust       = results.thrust_force_vector[0]

        # Calculate the aircraft drag for the clean configuration
        aerodynamics = SUAVE.Analyses.Aerodynamics.Fidelity_Zero()
        aerodynamics.geometry = vehicle
        aerodynamics.initialize()

        aero_clean = aerodynamics.evaluate(state)

        i = 0
        induced_drag_flap  = np.zeros(len(vehicle.wings.keys()))
        induced_drag       = np.zeros(len(vehicle.wings.keys()))
        parasite_drag_flap = np.zeros(len(vehicle.wings.keys()))
        interfer_drag_flap = np.zeros(len(vehicle.wings.keys()))
        parasite_drag_slat = np.zeros(len(vehicle.wings.keys()))
        interfer_drag_slat = np.zeros(len(vehicle.wings.keys()))
        AR_data            = np.zeros(len(vehicle.wings.keys()))
        component_lift     = np.zeros(len(vehicle.wings.keys()))
        #e                  = state.conditions.aerodynamics.drag_breakdown.induced.efficiency_factor[0,0,0,0]
        e                  = state.conditions.aerodynamics.drag_breakdown.induced.efficiency_factor[0]

        for wing in vehicle.wings:
            AR_data[i]        = wing.aspect_ratio
            component_lift[i] = aero_clean.lift.inviscid_wings[wing.tag][0,0]

            induced_drag_flap[i] = compute_induced_drag_highlift(wing,zero_aoa_lift_coef[i])
            induced_drag[i]      = (component_lift[i] + induced_drag_flap[i])**2/(np.pi*e*AR_data[i])
            parasite_drag_flap[i],parasite_drag_slat[i] = compute_parasite_drag_highlift(aero_clean.drag.parasite.wings[wing.tag], \
                                                          wing,alphad0[i],c_pr_c[i],flap_area[i], slat_area[i], Cla[i])
            interfer_drag_flap[i],interfer_drag_slat[i] = compute_interference_drag_highlift(wing,parasite_drag_flap[i],parasite_drag_slat[i])
            i += 1

        CDLG  = compute_landing_gear_drag(vehicle)
        CLtot = aero_clean.lift.total[0,0] + np.sum(zero_aoa_lift_coef)
        #CDtot = aero_clean.drag.total[0,0,0,0] - aero_clean.drag.induced[0,0,0,0] + phi*np.sum(induced_drag) + np.sum(parasite_drag_flap) +  \
        #        np.sum(parasite_drag_slat) + np.sum(interfer_drag_flap) + np.sum(interfer_drag_slat) + CDLG
        CDtot = aero_clean.drag.total[0] - aero_clean.drag.induced[0] + phi*np.sum(induced_drag) + np.sum(parasite_drag_flap) +  \
                np.sum(parasite_drag_slat) + np.sum(interfer_drag_flap) + np.sum(interfer_drag_slat) + CDLG

        Drag  = conditions.freestream.dynamic_pressure[0,0] * CDtot * reference_area
        Lift  = conditions.freestream.dynamic_pressure[0,0] * CLtot * reference_area

        accel     = (thrust[0] - Drag - surf_mu * (weight * sea_level_gravity - Lift - thrust[2]))/weight
        airspeed += accel * dtime
        
        ground_distance += airspeed_old*dtime + 0.5*accel*dtime**2
        # battery_draw    = conditions.propulsion.battery_energy[0,0]

        time_list.append(dtime)
        time_int_list.append(dtime)
        time += dtime

        '''conditions.freestream.dynamic_pressure  = 0.5 * rho * airspeed**2
        conditions.freestream.gravity[:,:]        = sea_level_gravity
        conditions.freestream.velocity[:,:]       = airspeed
        conditions.freestream.mach_number[:,:]    = airspeed/a
        conditions.freestream.speed_of_sound[:,:] = a
        conditions.freestream.temperature[:,:]    = T
        conditions.freestream.pressure[:,:]       = p
        conditions.freestream.reynolds_number[:,:] = rho*airspeed/mu'''
        conditions.freestream.dynamic_pressure = 0.5 * rho * airspeed**2
        conditions.freestream.gravity          = sea_level_gravity
        conditions.freestream.velocity         = airspeed
        conditions.freestream.mach_number      = airspeed/a
        conditions.freestream.speed_of_sound   = a
        conditions.freestream.temperature      = T
        conditions.freestream.pressure         = p
        conditions.freestream.reynolds_number  = rho*airspeed/mu

    # Distance to clear the obstacle
    rotation_distance   = 1.5 * V_LOF
    R_rotation          = 0.70842 * stall_speed**2 
    transition_distance = R_rotation * np.sin(np.arccos(1-h_obst/R_rotation))
    
    # Calculate takeoff distance
    takeoff_field_length = ground_distance + rotation_distance + transition_distance
    
    fid = open('SE2A_takeoff.dat','w')   # Open output file
    fid.write('Output file with Take-off details\n\n') #Start output printing
    fid.write( ' Take-off field length ................: ' + str( '%8.0F'   %   takeoff_field_length   ) + ' kg\n' )
    fid.close()
    
    # calculating second segment climb gradient, if required by user input
    #if compute_2nd_seg_climb:
        # Getting engine thrust at V2 (update only speed related conditions)
    '''state.conditions.freestream.dynamic_pressure  = np.array(np.atleast_1d(0.5 * rho * V2_speed**2))
    state.conditions.freestream.velocity          = np.array(np.atleast_1d(V2_speed))
    state.conditions.freestream.mach_number       = np.array(np.atleast_1d(V2_speed/ a))
    state.conditions.freestream.dynamic_viscosity = np.array(np.atleast_1d(mu))

    time_list     = [0, dtime]
    time_int_list = [dtime]
    state.numerics.time = time_list
    state.numerics.time = Data()
    state.numerics.time.integrate = time_int_list
    state.conditions.freestream.density           =  np.array(np.atleast_1d(rho))
    results = vehicle.propulsors.evaluate_thrust(state) # total thrust  
    thrust       = results.thrust_force_vector[0]

    # Compute windmilling drag
    windmilling_drag_coefficient = windmilling_drag(vehicle,state)

    # Compute asymmetry drag   
    asymmetry_drag_coefficient = asymmetry_drag(state, vehicle, windmilling_drag_coefficient)
           
    # Compute l over d ratio for takeoff condition, NO engine failure
    l_over_d = estimate_2ndseg_lift_drag_ratio(state,settings,vehicle) 
        
    # Compute L over D ratio for takeoff condition, WITH engine failure
    clv2 = maximum_lift_coefficient / (V2_VS_ratio) **2
    cdv2_all_engine = clv2 / l_over_d
    cdv2 = cdv2_all_engine + asymmetry_drag_coefficient + windmilling_drag_coefficient
    l_over_d_v2 = clv2 / cdv2
    
    # Compute 2nd segment climb gradient
    second_seg_climb_gradient = thrust / (weight*sea_level_gravity) - 1. / l_over_d_v2'''

    # return only takeoff_field_length
    return takeoff_field_length, battery_draw