# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 13:11:46 2020

@author: senth
"""

from scipy.interpolate import UnivariateSpline

def Auxillary_flap_eff(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Auxillary flap eff.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Flap_efficiency = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Flap_efficiency)
    
    return spl(flap_angle)

def Double_slotted_fixed_vane_eff(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Double slotted fixed vane eff.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Flap_efficiency = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Flap_efficiency)
    
    return spl(flap_angle)

def Triple_slotted_eff(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Triple slotted eff.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Flap_efficiency = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Flap_efficiency)
    
    return spl(flap_angle)

def Simple_hinge_eff(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Simple_hinge_eff.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Flap_efficiency = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Flap_efficiency)
    
    return spl(flap_angle)

def Fowler_eff(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fowler_eff.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Flap_efficiency = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Flap_efficiency)
    
    return spl(flap_angle)

def Fixed_hinge_1(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fixed_hinge_Cex.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Fixed_hinge_2(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fixed_hinge_Cex2.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Single_Slotted_Cex(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Single_Slotted.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Double_Slotted_fixed_vane_Cex(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Double_Slotted_fixed_vane.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Double_Slotted_variable_Cex(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Double_Slotted_variable.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Fowler_single_double_slotted_fixed_vane_Cex(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fowler_single_double_slotted_fixed_vane.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Fowler_double_triple_slotted_variable_Cex(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fowler_double_triple_slotted_variable.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Kb_(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Fowler_double_triple_slotted_variable.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Profiledragcf_c01(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Profiledragcf_c0.1.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Profiledragcf_c02(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Profiledragcf_c0.2.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Profiledragcf_c03(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Profiledragcf_c0.3.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def Profiledragcf_c04(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/Profiledragcf_c0.4.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def V_factor_A10_y100(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y1.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def V_factor_A10_y025(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y025.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def V_factor_A10_y050(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y050.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def V_factor_A10_y075(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y075.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def W_factor_A10_y100(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y1.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def W_factor_A10_y025(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y025.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def W_factor_A10_y050(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y050.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)

def W_factor_A10_y075(flap_angle):

    from numpy import genfromtxt
    my_data = genfromtxt('./curve_fit/V_factor_A10_y075.csv', delimiter=',')
    
    Flap_angle = my_data[:,0]
    Chord_extension_ratio = my_data[:,1]
    
    spl = UnivariateSpline(Flap_angle,Chord_extension_ratio)
    
    return spl(flap_angle)