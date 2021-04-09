## @ingroupMethods-Aeroelasticity-Fidelity_Zero
# NASA_Flutter.py
# 
# Created:  Jul 2020, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import math as m
from SUAVE.Core import Units, Data

# ----------------------------------------------------------------------
#  Atmospheric Attenuation
# ----------------------------------------------------------------------

## @ingroupMethods-Aeroelasticity-Fidelity_Zero
def calculate_flutter_speed(limit_criteria,vehicle):
    """ Calculate a flutter mach number for a conventional wing
            Inputs:
                wing.AR                                   [Unitless]
                     S                     - semi-span    [m]
                     taper                                [Unitless]
                     chords.root                          [m]
                     chords.tip                           [m]
                     sweeps.quarter_chord                 [deg]
                     materials.GJ_root                    
                     materials.GJ_mid 
                     materials.density                    [kg/cu m]
                weight_breakdown.wing                     [kg]


            Outputs: 
                M_flutter - Flutter Mach number           [Unitless]

            Assumptions:
                Present flutter analysis is applicable only for subsonic and transonic speeds"""
      

    # unpack
    MTOW          = vehicle.mass_properties.max_takeoff
    wing          = vehicle.wings['main_wing']
    AR            = 0.5 * wing.aspect_ratio
    b             = wing.spans.projected
    S             = 0.5 * wing.spans.projected
    sweep         = wing.sweeps.quarter_chord
    taper         = wing.taper
    root          = wing.chords.root
    GJ_root       = wing.aeroelasticity.GJ
    GJ_mid        = wing.aeroelasticity.GJ_mid
    CG            = wing.aeroelasticity.CG[0]
    d_cgr_i       = wing.aeroelasticity.CG_offset
    W_ex          = 0.5 * vehicle.weight_breakdown.wing
    x0_60         = wing.aeroelasticity.x0_60
    fspar_60      = wing.aeroelasticity.fspar_60
    rspar_60      = wing.aeroelasticity.rspar_60

    F_low_lim = limit_criteria.F_low_lim
    R_low_lim = limit_criteria.R_low_lim

    # Constants
    rho0 = 0.07657

    # Convert all values into US units
    MTOW = MTOW / Units.lb
    b    = b / Units.ft
    S    = S / Units.ft
    root = root / Units.ft
    GJ_root = GJ_root / 0.41348
    GJ_mid  = GJ_mid / 0.41348
    W_ex    = W_ex /Units.lb


    if len(wing.Segments.keys()) == 2:
        root          = root*wing.Segments[0].root_chord_percent
        taper         = wing.Segments[1].root_chord_percent
        sweep         = wing.Segments[0].sweeps.quarter_chord
        C_60          = 0.6*root*(1+taper)+root
    elif len(wing.Segments.keys()) > 2:
        print('The wing does not have straight leading and trailing edges or has more than 2 sections defined \n \
               Consider using a simple rectangular wing definition or a limit Segments to 2')

    C_60 = 0.60*root*(taper-1)+root
    C_75 = 0.75*root*(taper-1)+root

    # Check what high-lift divices exist and estimate the strip weight at 60% chord
    try:
        slats = wing.slats.type
    except NameError:
        slats = None

    if slats == 'hinged_nose' and wing.slats.span_end >= 0.6:
        W_slat = 0.204816 * (6.3935*np.log(MTOW)-32.674) * C_60 * fspar_60
    elif slats == 'slat' and wing.slats.span_end >= 0.6:
        W_slat = 0.204816 * (4.4001*np.log(MTOW)-23.927) * C_60 * fspar_60
    else:
        W_slat = 0

    try:
        flaps = wing.flaps.type
    except NameError:
        flaps = None

    if (flaps == 'single_slotted' or flaps == 'split') and wing.flaps.span_end >= 0.6:
        W_flap = 0.204816 * (4.3082E-5*MTOW+4.5865) * C_60 * (1-rspar_60)  
    elif flaps == 'double_slotted' and wing.flaps.span_end >= 0.6:
        W_flap = 0.204816 * (4.1142E-5*MTOW+12.837) * C_60 * (1-rspar_60)     
    elif flaps == 'Fowler' and wing.flaps.span_end >= 0.6:
        W_flap = 0.204816 * (2.8481E-5*MTOW+26.082) * C_60 * (1-rspar_60)  
    elif flaps == 'double_slotted_Fowler' and wing.flaps.span_end >= 0.6:
        W_flap = 0.204816 * 3.1884*MTOW**0.23242 * C_60 * (1-rspar_60)    
    elif flaps == 'triple_slotted_Fowler' and wing.flaps.span_end >= 0.6:
        W_flap = 0.204816 * 3.7225*MTOW**0.23013 * C_60 * (1-rspar_60)  
    else:
        W_flap = 0
    

    # Calculate running weight (assume that the running weight is proportional to the wing chord station)
    # Assume, the wing half CG is located at 35% span
    Ixx = W_ex * (CG*S)**2/S
    tip     = taper * root
    W2      = 2*W_ex/(S*(root/tip+1))
    Ixx2    = 2*Ixx/(S*(root/tip+1))
    W_60    = C_60/tip*W2                             # running weight at 60% semispan
    Ixx_60  = C_60/tip*Ixx2                     # running pitching moment of inertia
    Rgyb_60 = 2/C_60 * (Ixx_60/W_60)**0.5           # radius of gyration at 60% semispan

    # Calculate torsional frequency
    L      = S / np.cos(sweep) 
    Ka     = 2.9216e-1 * np.log(GJ_mid/GJ_root) + 1.5722  
    wa     = Ka/L * (GJ_root/(Ixx_60/32.2))**0.5   

    # Calculate the mass ratio
    miu0 = W_ex / (np.pi*rho0*(1+taper+taper**2)*root**2*S/12)
 
    # Comute CG location at 60% span wrt the chord
    Xflap = rspar_60 + 0.3 * (1-rspar_60) 
    Xslat = 0.5 * fspar_60 
    cgr_i = x0_60 - (W_flap*(x0_60 - Xflap) + W_slat*(x0_60 - Xslat))/(W_60 + W_flap + W_slat) + d_cgr_i

    # Compute correction factors for the actual wing
    K_AR    = 1 + 5/15 * (1/AR - 0.5)
    K_taper = -7.0789e-1 * taper**3 + 2.7689 * taper**2 - 3.2785 * taper + 2.1387
    K_miu   = -3.3505e-7 * miu0**3 + 9.8965e-5 * miu0**2 - 1.0268e-2 * miu0 + 1.232
    K_CG    = -7.0732e1 * cgr_i**3 + 1.2082e2 * cgr_i**2 - 6.9755e1 * cgr_i + 1.4439e1
    K_rgyb  = 1.3103 * Rgyb_60 + 0.3452

    K_all = K_AR * K_taper * K_miu * K_CG * K_rgyb

    # Compute the Regier number
    V_R = 0.5 * C_75 * wa * (miu0)**0.5
    R   = V_R/1116.45

    # Correct semi-empirical 
    F_low_lim = F_low_lim * K_all
    R_low_lim = R_low_lim / K_all

    # Interpolate the data to obtain the flutter Mach number
    F_lim     = np.interp(R,R_low_lim[:, 1],F_low_lim[:, 1])

    M_flutter = 0
    for i in range(np.shape(F_low_lim)[0]-1):
        if F_lim < F_low_lim[i, 1] and F_lim > F_low_lim[i+1, 1]:
            M_flutter = 0.5*(F_low_lim[i, 0] + F_low_lim[i+1, 0])

    if M_flutter == 0:
        print('Warning: the divergence Mach number is supersonic')         
    
    if sweep / Units.deg > 20:
        print('Warning: the sweep angle is larger than the model is built for\n \
               Use different surrogates') 

    # Calculate flutter equivalent airspeed
    V_flutter = M_flutter * V_R

    # Pack results
    flutter = Data()
    flutter.M = M_flutter
    flutter.V = V_flutter

    return flutter
