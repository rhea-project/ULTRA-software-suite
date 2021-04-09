# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:55:35 2020

@author: senth
"""

import math


def Parasite_drag(vehicle,Mu,air_temperature,reference_airspeed, air_density, transaction_point,IF):
    
   
    
    Mu = (1.458*10**-6)*(air_temperature**1.5)*(1/(air_temperature+110.4))

    
    reference_area = vehicle.reference_area
    
    
    #Parasite main wing
    
    refrence_length = vehicle.wings.main_wing.chords.root
    
    Re = air_density*refrence_length*reference_airspeed/Mu
    
    K = 2.088*10**-5
    
    Re2 = 38.21*(refrence_length/K)**1.053
    
    if Re>Re2:
        Re = Re2
    
    airfoil_thickness_ratio = vehicle.wings.main_wing.thickness_to_chord
    Sweepangle_midchord = vehicle.wings.main_wing.sweeps.quarter_chord
    
    FFwing = (3.3*airfoil_thickness_ratio - 0.008*airfoil_thickness_ratio**2 + 27.0*airfoil_thickness_ratio**3)*math.cos(Sweepangle_midchord)**2 +1
    
    Xtr = transaction_point[1]
    X0 = transaction_point[2]
    
    Cf = (0.074/(Re**0.2))*(1-((Xtr-X0)/refrence_length))**0.8
    
    Splanform = vehicle.wings.main_wing.areas.reference
    
    SWet = 2.05*Splanform
    
    CdWing = (1/reference_area)*Cf*FFwing*IF*SWet
    
    #Parasite Horizontal Stabilizer
    
    refrence_length = vehicle.wings.horizontal_stabilizer.chords.root
    
    Re = air_density*refrence_length*reference_airspeed/Mu
    
    K = 2.088*10**-5
    
    Re2 = 38.21*(refrence_length/K)**1.053
    
    if Re>Re2:
        Re = Re2
    
    airfoil_thickness_ratio = vehicle.wings.horizontal_stabilizer.thickness_to_chord
    Sweepangle_midchord = vehicle.wings.horizontal_stabilizer.sweeps.quarter_chord
    
    FFtail = (3.52*(airfoil_thickness_ratio))*math.cos(Sweepangle_midchord)**2 +1
    
    Xtr = transaction_point[3]
    X0 = transaction_point[4]
    
    Cf = (0.074/(Re**0.2))*(1-((Xtr-X0)/refrence_length))**0.8
    
    Splanform = vehicle.wings.horizontal_stabilizer.areas.reference
    
    SWet = 2.05*Splanform
    
    Cdtail = (1/reference_area)*Cf*FFtail*IF*SWet
    
    #Parasite Fuselage
        
    refrence_length = vehicle.fuselages.fuselage.lengths.total
    
    Re = air_density*refrence_length*reference_airspeed/Mu
    
    K = 2.088*10**-5
    
    Re2 = 38.21*(refrence_length/K)**1.053
    
    if Re>Re2:
        Re = Re2
    
    Xtr = transaction_point[5]
    X0 = transaction_point[6]
    
    Cf = (0.074/(Re**0.2))*(1-((Xtr-X0)/refrence_length))**0.8
    
    Amax = vehicle.fuselages.fuselage.areas.front_projected
    
    length_of_body = vehicle.fuselages.fuselage.lengths.total
    
    f = length_of_body/math.sqrt(4*Amax/math.pi)
    
    FFbody = 1 + (2.2/f**1.5) + (3.8/f**3)
    
    Cf = (0.074/(Re**0.2))*(1-((Xtr-X0)/refrence_length))**0.8
    
    SWet = vehicle.fuselages.fuselage.areas.wetted
    
    Cdbody = (1/reference_area)*Cf*FFbody*IF*SWet
    
    Cdmin = CdWing + Cdtail + Cdbody
    
    return Cdmin