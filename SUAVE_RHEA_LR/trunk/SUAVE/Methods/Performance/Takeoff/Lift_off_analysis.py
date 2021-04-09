# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:37:57 2020

@author: senth
"""

# 
# 
# Created:  
# Modified: 

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# Python Imports
import numpy as np
import pylab as plt
import math
#import matplotlib.pyplot as plt

# SUAVE Imports

import SUAVE
from SUAVE.Core import Data, Units

from SUAVE.Attributes.Atmospheres.Atmosphere import Atmosphere
from SUAVE.Analyses.Mission.Segments.Conditions import Aerodynamics,Numerics

from Parasite_drag import Parasite_drag
from  Lift_induced_drag import  Lift_induced_drag
from miscellaneous_drag_aircraft import miscellaneous_drag_aircraft
from Flap_lift import Flap_lift
# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def Lift_off_analysis(vehicle, conditions, airfoil_lift):
    
    T_max           = 40* Units.s
    del_T           = 1 * Units.s
    
    AOA = 2 * Units.deg
    AOA0 = 0 * Units.deg 
    
    # ==============================================
    # Computing atmospheric conditions
    # ==============================================
    planet = SUAVE.Analyses.Planets.Planet()
    atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    
    altitude        = 0 * Units.ft
    delta_isa       = 0
    
    atmo_values       = atmosphere.compute_values(altitude,delta_isa)
    conditions        = SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics()
    
    p   = atmo_values.pressure
    T   = atmo_values.temperature
    rho = atmo_values.density
    a   = atmo_values.speed_of_sound
    mu  = atmo_values.dynamic_viscosity
    sea_level_gravity = atmosphere.planet.sea_level_gravity 
    #reynolds_number = atmo_values.reynolds_number 
    
    # ==============================================
    # Unpack
    # ==============================================
    
    Thrust          = vehicle.propulsors.turbofan.thrust.total_design
    Weight          = vehicle.mass_properties.takeoff*sea_level_gravity
    Mass            = vehicle.mass_properties.takeoff 
    wing_area       = vehicle.wings.main_wing.areas.reference
    
    
    try:
        V2_VS_ratio = vehicle.V2_VS_ratio
    except:
        V2_VS_ratio = 1.20
    
    # ==============================================
    # Determining vehicle maximum lift coefficient
    # ==============================================
    try:   # aircraft maximum lift informed by user
        maximum_lift_coefficient = vehicle.maximum_lift_coefficient
    except:
        # Using semi-empirical method for maximum lift coefficient calculation
        from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Lift import compute_max_lift_coeff

        # Condition to CLmax calculation: 90KTAS @ 10000ft, ISA
        conditions  = atmosphere.compute_values(10000. * Units.ft)
        conditions.freestream=Data()
        conditions.freestream.density   = conditions.density
        conditions.freestream.dynamic_viscosity = conditions.dynamic_viscosity
        conditions.freestream.velocity  = 90. * Units.knots
        try:
            maximum_lift_coefficient, induced_drag_high_lift = compute_max_lift_coeff(vehicle,conditions)
            vehicle.maximum_lift_coefficient = maximum_lift_coefficient
        except:
            raise ValueError("Maximum lift coefficient calculation error. Please, check inputs")
    
    
     # ==============================================
    # Computing speeds (Vs, V2, 0.7*V2)
    # ==============================================
    stall_speed = (2 * Weight * sea_level_gravity / (rho * wing_area * maximum_lift_coefficient)) ** 0.5
    V2_speed    = V2_VS_ratio * stall_speed
    speed_for_thrust  = 0.70 * V2_speed

    # ==============================================
    # Determining vehicle number of engines
    # ==============================================
    engine_number = 0.
    for propulsor in vehicle.propulsors : # may have than one propulsor
        engine_number += propulsor.number_of_engines
    if engine_number == 0:
        raise ValueError("No engine found in the vehicle")

    # ==============================================
    # Getting engine thrust
    # ==============================================    
    state = Data()
    state.conditions = Aerodynamics() 
    state.numerics   = Numerics()
    conditions = state.conditions    

    conditions.freestream.dynamic_pressure = np.array(np.atleast_1d(0.5 * rho * speed_for_thrust**2))
    conditions.freestream.gravity          = np.array([np.atleast_1d(sea_level_gravity)])
    conditions.freestream.velocity         = np.array(np.atleast_1d(speed_for_thrust))
    conditions.freestream.mach_number      = np.array(np.atleast_1d(speed_for_thrust/ a))
    conditions.freestream.speed_of_sound   = np.array(a)
    conditions.freestream.temperature      = np.array(np.atleast_1d(T))
    conditions.freestream.pressure         = np.array(np.atleast_1d(p))
    conditions.propulsion.throttle         = np.array(np.atleast_1d(1.))
    conditions.freestream.reynolds_number  = np.array(np.atleast_1d(12e6))

    
    results = vehicle.propulsors.evaluate_thrust(state) # total thrust
    
    thrust = np.array(results.thrust_force_vector)
    Thrust = thrust[0][0]
    #print(Thrust)
       
    #Thrust = 48000
    # ==============================================
    
    AR = vehicle.wings.main_wing.aspect_ratio
    beta = 1
    AC2 = vehicle.wings.main_wing.sweeps.quarter_chord 
    CL_alpha = .1
    kappa= CL_alpha*(180/math.pi)/(2*math.pi)

    # ==============================================
    # calculate lift coefficient
    # ==============================================
    gamma = 15 * Units.deg
    Lift_Coefficient   = Calculate_Lift_Coefficient(AR, beta, kappa, AC2, AOA, AOA0)
    LCa_gamma = Calculate_Lift_Coefficient(AR, beta, kappa, AC2, gamma, AOA0)
    Flap_lift_coefficient = Flap_lift(vehicle,gamma,Lift_Coefficient,LCa_gamma, airfoil_lift)
    
    Total_lift_coefficient = Lift_Coefficient + Flap_lift_coefficient
    
    N = T_max/del_T +1
    N = int(N)
    X_Velocity = np.zeros(N) 
    Y_Velocity = np.zeros(N) 
    X_Distance = np.zeros(N) 
    Y_Distance = np.zeros(N)  
    
    transaction_point = np.zeros(7)
    
    transaction_point[1] = 4 * Units.meter
    transaction_point[2] = 3 * Units.meter
    
    transaction_point[3] = 2.5 * Units.meter
    transaction_point[4] = 2 * Units.meter
    
    transaction_point[5] = 20 * Units.meter
    transaction_point[6] = 15 * Units.meter
    
    IF = 1.0
    


    Cdmin = Parasite_drag(vehicle,mu, conditions.freestream.temperature ,conditions.freestream.velocity, rho, transaction_point,IF)
    
    Oswald_efficiency = 0.9
    Cdi = Lift_induced_drag(Total_lift_coefficient, AR , Oswald_efficiency)
    
    Cdm = miscellaneous_drag_aircraft(conditions,vehicle)
    
    
    

    
    Drag_Coefficient   = Cdmin + Cdi + Cdm
    
    n=0
    for T in range(0,T_max,del_T):
        
        Cdmin = Parasite_drag(vehicle,mu, conditions.freestream.temperature ,conditions.freestream.velocity + X_Velocity[n] , rho, transaction_point,IF)
    
        Oswald_efficiency = 0.9
        
        Cdi = Lift_induced_drag(Lift_Coefficient, AR , Oswald_efficiency)
        
        Drag_Coefficient   = Cdmin + Cdi + Cdm
        
        Drag = Lift_or_Drag(X_Velocity[n], rho, wing_area , Drag_Coefficient)
        
        X_Force = Thrust - Drag
        
        X_Acceleration = X_Force/Mass
        
        X_Velocity[n+1] = X_Velocity[n] + (X_Acceleration*del_T**2)/2
        
        X_Distance[n+1] = X_Distance[n] + X_Velocity[n]*del_T + (1/2)*X_Acceleration*del_T**2
        
        M = X_Velocity[n]/a
        
        #print(X_Velocity[n]*3.6)
        if (n > 30):
            AOA = 20 * Units.deg
            
        beta = math.sqrt(1-M**2)
        
        kappa= Lift_Coefficient*(180/math.pi)/(2*math.pi)
        
        Lift_Coefficient = Calculate_Lift_Coefficient(AR, beta, kappa, AC2, AOA, AOA0)
        
        LCa_gamma = Calculate_Lift_Coefficient(AR, beta, kappa, AC2, gamma, AOA0)
        Flap_lift_coefficient = Flap_lift(vehicle,gamma,Lift_Coefficient,LCa_gamma, airfoil_lift)      
        Total_lift_coefficient = Lift_Coefficient + Flap_lift_coefficient
        #print(Total_lift_coefficient)
        Lift = Lift_or_Drag(X_Velocity[n], rho, wing_area , Total_lift_coefficient)
        
        Y_Force = Lift - Weight
        
        if (Y_Force < 0):
            Y_Force = 0
        
        Y_Acceleration = Y_Force/Mass
        
        Y_Velocity[n+1] = Y_Velocity[n] + (Y_Acceleration*del_T**2)/2
        
        Y_Distance[n+1] = Y_Distance[n] + Y_Velocity[n]*del_T + (1/2)*Y_Acceleration*del_T**2
        
        n=n+1
        
    
    plt.figure(1)
    plt.plot(X_Distance,Y_Distance)


def Lift_or_Drag(velocity, rho, S , Lift_Coefficient_or_Drag_Coefficient):
    Lift_or_Drag = rho*(velocity**2)*S*Lift_Coefficient_or_Drag_Coefficient/2
    
    return Lift_or_Drag


def Calculate_Lift_Coefficient(AR, beta, kappa, AC2, AOA, AOA0):
    Lift_Coefficient = (2*math.pi*AR)/(2+ math.sqrt( ( (AR*beta/kappa)**2 * ( 1+((math.tan(AC2)**2)/beta**2)) ) +4 ))    
    Lift_Coefficient = Lift_Coefficient*(AOA-AOA0)
    return Lift_Coefficient

