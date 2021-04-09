## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Drag
# induced_drag_aircraft.py
# 
# Created:  Dec 2013, SUAVE Team
# Modified: Jan 2016, E. Botero
#                     S. Karpuk
       

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# suave imports
from SUAVE.Core import Data, Units

# package imports
import math  as m
import numpy as np

# ----------------------------------------------------------------------
#  Induced Drag Aircraft
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Drag
def induced_drag_aircraft(state,settings,geometry):
    """Determines induced drag for the full aircraft

    Assumptions:
    Based on fits

    Source:
    adg.stanford.edu (Stanford AA241 A/B Course Notes)
    M. Nita, D. Scholz, 'Estimating the Oswald Factor from Basic Aircraft Geometrical Parameters'
        hamburg University of Applied Sciences, Aero - Aircraft Design and Systems Group

    Inputs:
    state.conditions.aerodynamics.lift_coefficient               [Unitless]
    state.conditions.aerodynamics.drag_breakdown.parasite.total  [Unitless]
    configuration.oswald_efficiency_factor                       [Unitless]
    configuration.viscous_lift_dependent_drag_factor             [Unitless]
    geometry.wings['main_wing'].span_efficiency                  [Unitless]
    geometry.wings['main_wing'].aspect_ratio                     [Unitless]
    geometry.wings['main_wing'].spans.projected                  [m]
    geometry.wings['main_wing'].sweeps.quarter_chord             [rad]
    geometry.wings['main_wing'].taper                            [Unitless]
    geometry.fuselages['fuselage'].width                         [m]

    Outputs:
    total_induced_drag                                           [Unitless]

    Properties Used:
    N/A
    """

    # unpack inputs
    conditions    = state.conditions
    configuration = settings
    
    aircraft_lift = conditions.aerodynamics.lift_coefficient
    mach          = conditions.freestream.mach_number
    e             = configuration.oswald_efficiency_factor
    K             = configuration.viscous_lift_dependent_drag_factor
    span          = geometry.wings['main_wing'].spans.projected 
    ar            = geometry.wings['main_wing'].aspect_ratio 
    CDp           = state.conditions.aerodynamics.drag_breakdown.parasite.total
    CDc           = state.conditions.aerodynamics.drag_breakdown.compressible.total
    t_c_w         = geometry.wings['main_wing'].thickness_to_chord
    taper         = geometry.wings['main_wing'].taper
    sweep         = geometry.wings['main_wing'].sweeps.quarter_chord / Units.degrees
    sweep_deg     = sweep * Units.degrees
    
    if 'fuselage' in geometry.fuselages:
        d_f = geometry.fuselages['fuselage'].width
    else:
        d_f = 0

    cl_w      = 0  
    cos_sweep = np.cos(sweep_deg)

    # get effective Cl and sweep
    tc = t_c_w /(cos_sweep)
    cl = cl_w / (cos_sweep*cos_sweep)

    # compressibility drag based on regressed fits from AA241
    mcc_cos_ws = 0.922321524499352       \
               - 1.153885166170620*tc    \
               - 0.304541067183461*cl    \
               + 0.332881324404729*tc*tc \
               + 0.467317361111105*tc*cl \
               + 0.087490431201549*cl*cl
        
    # crest-critical mach number, corrected for wing sweep
    mcc = mcc_cos_ws / cos_sweep

    # divergence ratio
    mo_mc = mach/mcc
    
    # compressibility correlation, Shevell
    dcdc_cos3g = 0.0019*mo_mc**14.641
    
    # compressibility drag
    CDc = dcdc_cos3g * cos_sweep*cos_sweep*cos_sweep

    if e == None:
        dtaper = -0.357+0.45*m.exp(-0.0375*np.abs(sweep)) 
        s      = 1 - 2 * (d_f/span)**2
        f      = 0.0524*(taper-dtaper)**4 - 0.15*(taper-dtaper)**3 + \
                 0.1659*(taper-dtaper)**2 - 0.0706*(taper-dtaper) + 0.0119
        u      = 1/(1+f*ar) 
        e      = 1/(1/(u*s)+np.pi*ar*K*(CDp+CDc))

    # start the result
    total_induced_drag = aircraft_lift**2 / (np.pi*ar*e)

    # For SBW case, reduce the strut's CDi， -0.0008 comes from SUGAR case
    total_induced_drag = total_induced_drag - 0.0008
        
    # store data
    conditions.aerodynamics.drag_breakdown.induced = Data(
        total             = total_induced_drag ,
        efficiency_factor = e                  ,
        aspect_ratio      = ar                 ,
    )

    return total_induced_drag
