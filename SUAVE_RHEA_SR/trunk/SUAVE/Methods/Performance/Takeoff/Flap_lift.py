# -*- coding: utf-8 -*-

from Curve_fit import *
import math
from SUAVE.Core import Data, Units

def Flap_lift(vehicle,flap_angle,CLa,LCa_gamma, airfoil_lift):
    
    flap_angle = 30 * Units.deg
    
    if vehicle.wings.main_wing.flaps.type == 'double_slotted':
        eff = Double_slotted_fixed_vane_eff(flap_angle)
        Cex = Double_Slotted_fixed_vane_Cex(flap_angle)
        incr = 2
        
    elif vehicle.wings.main_wing.flaps.type == 'single_slotted':
        eff = Simple_hinge_eff(flap_angle)
        Cex = Single_Slotted_Cex(flap_angle)
        incr = 1
        
    elif vehicle.wings.main_wing.flaps.type == 'split':
        eff = Simple_hinge_eff(flap_angle)
        Cex = Fixed_hinge_1(flap_angle)
        Cex = Fixed_hinge_2(flap_angle)
        incr = 1
    
    elif vehicle.wings.main_wing.flaps.type == 'triple_slotted':
        eff = Triple_slotted_eff(flap_angle)
        Cex = Fowler_double_triple_slotted_variable_Cex(flap_angle)
        incr = 3
        
    elif vehicle.wings.main_wing.flaps.type == 'fowler':
        eff = Fowler_eff(flap_angle)
        Cex = Fowler_double_triple_slotted_variable_Cex(flap_angle)
        incr = 3
        
    else:
        print("invalid flaps type")
        
    Cf = vehicle.wings.main_wing.flaps.chord * vehicle.wings.main_wing.chords.mean_aerodynamic
    C  = vehicle.wings.main_wing.chords.mean_aerodynamic
    
    C_extended = C*(1+Cex*(Cf/C))
    
    Flap_span = (vehicle.wings.main_wing.flaps.span_end - vehicle.wings.main_wing.flaps.span_start) * vehicle.wings.main_wing.spans.projected
    
    theta = math.acos((2*Cf/C_extended)-1)
    alphaF = 1 - ((theta - math.sin(theta))/math.pi)
    
    del_Cl = eff*CLa*alphaF*flap_angle*incr
    Kb = 0.5
    Flap_increment = del_Cl*LCa_gamma*Kb/(airfoil_lift*eff*alphaF)
    
    return Flap_increment
    
        
    
        
    