## @ingroup Methods-Aerodynamics-Fidelity_Zero-Drag
# compute_landing_gear_drag.py
#
# Created:  Aug 2020, S. Karpuk
# Modified: 
#            

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Units
import numpy as np

# ----------------------------------------------------------------------
#  compute_landing_gear_drag
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_landing_gear_drag(vehicle):
    """Computes landing gear drag increment 

    Assumptions:
        Roskam's method for retractable landing gear only
    Source:
    Unknown

    Inputs:


    Outputs:


    Properties Used:
    N/A
    """          

    # Unpack inputs
    Sref      = vehicle.reference_area / Units['ft**2']

    Dmain        = vehicle.landing_gear.main_tire_diameter / Units.feet
    Dnose        = vehicle.landing_gear.nose_tire_diameter / Units.feet   
    Lstr_main    = vehicle.landing_gear.main_strut_length / Units.feet
    Lstr_nose    = vehicle.landing_gear.nose_strut_length / Units.feet
    N_main       = vehicle.landing_gear.main_units                                    # number of main landing gear units
    N_nose       = vehicle.landing_gear.nose_units                                    # number of nose landing gear
    Wmain        = vehicle.landing_gear.main_tire_width / Units.feet 
    Wnose        = vehicle.landing_gear.nose_tire_width / Units.feet 
    Sa_main      = vehicle.landing_gear.nose_gear_front_area_ratio 
    Sa_nose      = vehicle.landing_gear.nose_gear_front_area_ratio 

    Hmain = Lstr_main + 0.5*Dmain
    Hnose = Lstr_nose + 0.5*Dnose 
    if vehicle.landing_gear.closed is True:
        LG_coef = 0.04955
    else:
        LG_coef = 0.05328
    
    dCD_LG = N_nose * (LG_coef * np.exp(5.615*Sa_nose) * Hnose*Wnose/Sref) + \
             N_main * (LG_coef * np.exp(5.615*Sa_main) * Hmain*Wmain/Sref)
           
  
    return dCD_LG
