## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import math
from SUAVE.Core import Units
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Fuselage Weight Calculation
# ----------------------------------------------------------------------

def weight_fuselage(fues_length,cabin_length,wing_sweep,n_engines,CARGF,n_fues,max_width,
                    max_depth,DG,aircraft_type,cab_area,var_sw_fac,ULF,RSPSOB,p_h,p_0,VCMN):
    """ Calculate the  fuselage weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        fues_length  - Total fuselage length                 [m]
        cabin_length - Passenger cabin length                [m]
        wing_sweep   - Wing leading-edge sweep               [deg]
        n_engines    - Number of fuselage mounted engines                        
        CARGF        - Cargo aircraft floor factor:  0.0  for passenger transport
                                                     1.0  for military cargo transport
       
        aircraft_type                                        [Unitless]

        n_fues        - Number of fuselages      
        max_width     - Maximum fuselage width                [m]
        max_depth     - Maximum fuselage depth                [m]
        DG            - Design gross weight                   [kg]
        cab_area      - Passenger cabin floor area            [m²]
        var_sw_fac    - Wing variable sweep weight penalty factor: 0.0 for fixed-geometry wing 
                                                                   1.0 for full variable-sweep wing 
        ULF          - Structural ultimate load factor: default value 3.75      
        RSPSOB       - Percent chord of the HWB fuselage rear spar at the side of body: default =RSPCHD=70%        
        p_h          - Pressure at altitude H                 [Pa]
        p_0          - Pressure at sea level                  [Pa]
        VCMN         - Cruise Mach number
        
    Outputs:
      weight_fuselage                                        [kg]          
        
    Properties Used:
        N/A
    """ 
    
    # unpack inputs
    fues_length  = fues_length/Units.ft          # Convert meters to ft
    cabin_length = cabin_length/Units.ft         # Convert meters to ft
    max_width    = max_width/Units.ft            # Convert meters to ft 
    max_depth    = max_depth/Units.ft            # Convert meters to ft 
    DG           = DG/Units.lb                   # Convert kg to lb
    cab_area     = cab_area/Units.ft**2          # Convert meters squared to ft squared

    fix_side_wall_len = cabin_length - 0.5 * math.tan(wing_sweep) * max_width
    
    #For transport aircraft
    if aircraft_type =='transport':
        aver_diam             = (max_width+max_depth)/2
        w_fueselage_transport = 1.35*(fues_length*aver_diam)**1.28*(1+0.05*n_engines)*(1+0.38*CARGF)* n_fues
        weight_fueselage      = w_fueselage_transport
    #For HWB aircraft
    elif aircraft_type =='HWB':
        w_fueselage_HWB  = 1.8*DG**0.167*cab_area**1.06
        weight_fueselage = w_fueselage_HWB
    #For fighter/attack aircraft
    elif aircraft_type =='fighter':
        w_fueselage_fighter = 0.15*fues_length**0.9*DG**0.61*(1+0.3*n_fues)*(1+0.33*var_sw_fac)*n_fues**0.3
        weight_fueselage    = w_fueselage_fighter
    #For general aviation aircraft  
    else:
        delta              = p_h/p_0
        QCRUS              = 1481.35*delta*VCMN**2
        aver_diam          = (max_width+max_depth)/2
        fues_wetted_area   = 3.14159*((fues_length/aver_diam)-1.7)*aver_diam**2
        w_fueselage_gen_av = 0.052*fues_wetted_area**1.086*(ULF*DG)**0.177*QCRUS**0.241
        weight_fueselage   = w_fueselage_gen_av
        
    weight_fueselage = weight_fueselage*Units.lb
    
    return weight_fueselage

