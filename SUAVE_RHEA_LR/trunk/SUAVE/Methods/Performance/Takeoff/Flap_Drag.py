# -*- coding: utf-8 -*-
"""
Created on Tue Apr  25 15:59:23 2020

@author: senth
"""

from Curve_fit import *
import math

def Flap_profile_Drag(vehicle,flap_angle,CLa,CLo,Flap_lift_Coefficient, Drag_coefficient):
    
    Cf_by_C = vehicle.wings.main_wing.flaps.chord
    
    if(Cf_by_C <= 0.1):
        F = Profiledragcf_c01(flap_angle)
    elif((Cf_by_C <= 0.2) and (Cf_by_C > 0.1)):
        F = Profiledragcf_c02(flap_angle)
    elif((Cf_by_C <= 0.3) and (Cf_by_C > 0.2)):
        F = Profiledragcf_c03(flap_angle)
    elif((Cf_by_C <= 0.4) and (Cf_by_C > 0.3)):
        F = Profiledragcf_c02(flap_angle)
        
    Thickness_to_chord = vehicle.wings.main_wing.thickness_to_chord
    
    del_f_C_dpO = 0.55*Cf_by_C*F*pow((Cf_by_C/pow(Thickness_to_chord,(3/2))),(2/9))
    
    K2_3 = 1.15
    Kl   = 0.8
    
    Cf = vehicle.wings.main_wing.flaps.chord * vehicle.wings.main_wing.chords.mean_aerodynamic
    Flap_span = (vehicle.wings.main_wing.flaps.span_end - vehicle.wings.main_wing.flaps.span_start) * vehicle.wings.main_wing.spans.projected
    
    Area_Flap = Cf*Flap_span
    wing_area       = vehicle.wings.main_wing.areas.reference
    
    Swf_S = Area_Flap/wing_area
    
    AC2 = vehicle.wings.main_wing.sweeps.quarter_chord
    
    del_f_C_dp  = K2_3*Swf_S*del_f_C_dpO*math.cos(AC2)  -  Kl*Drag_coefficient*Flap_lift_Coefficient*(CLa -(CLo+Flap_lift_Coefficient/4) )
    
    return del_f_C_dp

def Flap_Vortex_Drag(vehicle,flap_angle,CLa,CLo,Flap_lift_Coefficient, Drag_coefficient):
    
    Kff = 1/3
    
    chord_tip = vehicle.wings.main_wing.chords.tip
    chord_root = vehicle.wings.main_wing.chords.root
    
    y = chord_tip/chord_root
    
    bf_b = vehicle.wings.main_wing.flaps.span_end - vehicle.wings.main_wing.flaps.span_start
    
    z = (0.07/(1+y))*pow((1-Kff),2)*bf_b
    
    if(bf_b <= 0.25):
        w = W_factor_A10_y025(bf_b)
        v = V_factor_A10_y025(bf_b)
    elif((bf_b > 0.25) and (bf_b <= 0.50)):
        w = W_factor_A10_y050(bf_b)
        v = V_factor_A10_y050(bf_b)
    elif((bf_b > 0.50) and (bf_b <= 0.75)):
        w = W_factor_A10_y075(bf_b)
        v = V_factor_A10_y075(bf_b)
    elif((bf_b > 0.50) and (bf_b <= 1)):
        w = W_factor_A10_y100(bf_b)
        v = V_factor_A10_y100(bf_b)
        
    del_f_CDv = (w+z)*pow(Flap_lift_Coefficient,2) + v*CLa*Flap_lift_Coefficient
    
    return del_f_CDv
    
    
    
    