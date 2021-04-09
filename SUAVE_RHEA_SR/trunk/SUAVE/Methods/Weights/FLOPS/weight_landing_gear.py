## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import math as m
import pandas as pd
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Landing Gear Weight Calculation
# ----------------------------------------------------------------------

def weight_lg(ac_type,engine_pos,wing_dih_angle,out_engine_pos,ac_landing_w,averg_diam_engine, \
              max_fues_width,tot_fues_length,design_range,len_ext_nose_lg, \
              len_ext_main_lg,cruise_type,ramp_w,CARBAS=0):
    
    """ Calculate the landing grear weight based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type         : aircraft type                             [Unitless]
        engine_pos      : engines position: 'on wings' or 'otherwise'
                                                                             
        wing_dih_angle  : Wing dihedral angle                       [deg]  
        out_engine_pos  : Location of outboard engine               [m] 
        ac_landing_w    : Aircraft design landing weight            [kg]
        len_ext_main_lg : Length of the extended main landing gear  [m]
        CARBAS          : Carrier based aircraft switch: 1.0 is for carrier-based aircraft  
                                                        0.0 is for land-based aircraft
        len_ext_nose_lg : Length of the extended nose landing gear  [m]
        ramp_w          : Ramp weight                               [kg] #not needed if the landing weight is given
        cruise_type     : 1 for subsonic
                          2 for supersonic
        design_range                                                [nmi] #not needed if the landing weight is given
        averg_diam_engine : Average diameter of each scaled engine, scaled to account for distributed propulsion if applicable [ft]
        max_fues_width  :  Maximum fuselage width                   [m]
        tot_fues_length :  Total fuselage length                    [m]  #not needed if the engine is wing-mounted
    Outputs:
       total_lg_weight                                              [kg]
                                                 
        
    Properties Used:
        N/A
    """ 
    # unpack inputs
    design_range      = 0.000539957*design_range
    out_engine_pos    = out_engine_pos/Units.inch   # Convert meter to inches
    ac_landing_w      = ac_landing_w/Units.lb       # Convert kg to lb
    max_fues_width    = max_fues_width/Units.ft
    tot_fues_length   = tot_fues_length/Units.ft
    averg_diam_engine = averg_diam_engine/Units.ft
    len_ext_nose_lg   = len_ext_nose_lg/Units.inch
    len_ext_main_lg   = len_ext_main_lg/Units.inch
    ramp_w            = ramp_w/Units.lb

    #Calculate the landing gear weight using FLOPS methods
    if cruise_type==1: 
        cruise_factor = 0.00004
    elif  cruise_type==2:
        cruise_factor = 0.00009
    #Calculate the length of the extended main landing gear oleo
    if len_ext_main_lg==0:
        if engine_pos=='on wings':
            len_ext_main_lg = 12*averg_diam_engine+(0.26-m.tan(wing_dih_angle))*(out_engine_pos-6*max_fues_width) 
        else:
             len_ext_main_lg = 0.75*tot_fues_length

      #Calculate the length of the extended nose landing gear oleo
    if len_ext_nose_lg==0:
       len_ext_nose_lg = 0.7*len_ext_main_lg
    if ac_landing_w==0:
        ac_landing_w = ramp_w*(1-cruise_factor*design_range)

        if ac_type == 'fighter':
           ac_type_ind = 1
        else:
            ac_type_ind = 0
            
    #Calculate the weight of the nose landing gear
    # 0.85 is for future aircraft, comes from SUGAR_REPORT_2011
    weight_nose_landing_gear = (0.048-0.008*ac_type_ind)*ac_landing_w**0.67* len_ext_nose_lg**0.43*(1+0.8*CARBAS) * 0.85
    #Calculate the weight of the main landing gear
    # 0.85 is for future aircraft, comes from SUGAR_REPORT_2011
    weight_main_landing_gear = (0.0117-0.0012*ac_type_ind)*ac_landing_w**0.95*len_ext_main_lg**0.43 * 0.85   
        
    total_lg_weight          = weight_main_landing_gear+weight_nose_landing_gear
    total_lg_weight          = total_lg_weight*Units.lb
    weight_nose_landing_gear = weight_nose_landing_gear*Units.lb
    weight_main_landing_gear = weight_main_landing_gear*Units.lb

    return total_lg_weight, weight_nose_landing_gear, weight_main_landing_gear  


