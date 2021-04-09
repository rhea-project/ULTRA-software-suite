# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:55:35 2020

@author: senth
"""

import math

def Skin_friction():
    
    air_temperature 
    
    air_density
    
    refrence_length
    
    reference_airspeed
    
    FF
    IF
    
    Xtr
    X0
    
    
    Mu = (1.458*10**-6)*(air_tempareture**1.5)(1/(air_temperature+110.4))
    
    Re = density*refrence_length*reference_airspeed/Mu
    
    K = 2.08810**-5
    
    Re2 = 38.21(refrence_length/K)**1.053
    
    if Re>Re2:
        Re = Re2
    
    airfoil_thickness_ratio
    Sweepangle_midchord
    
    FFwing = (3.3*airfoil_thickness_ratio - 0.008*airfoil_thickness_ratio**2 + 27.0*airfoil_thickness_ratio**3)*cos(Sweepangle_midchord)**2 +1
    
    FFtail = (3.52*(airfoil_thickness_ratio))*cos(Sweepangle_midchord)**2 +1
    
    Amax
    length_of_body
    f = length_of_body/math.sqrt(4*Amax/math.pi)
    
    FFbody = 1 + (2.2/f**1.5) + (3.8/f**3)
    
    Cf = (0.074/(Re**0.2))(1-((Xtr-X0)/refrence_length))**0.8
    
    SWet = 2.05*Splanform
    
    Cd = (1/refrence_length)*Cf*FF*IF*SWet
    
    return Cd