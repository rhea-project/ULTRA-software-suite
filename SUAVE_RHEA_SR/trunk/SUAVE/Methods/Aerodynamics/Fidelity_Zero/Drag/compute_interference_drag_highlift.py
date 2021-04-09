## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
# compute_flap_lift.py
#
# Created:  Dec 2013, A. Varyar
# Modified: Feb 2014, T. Orra
#           Jan 2016, E. Botero         

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Units
import numpy as np

# ----------------------------------------------------------------------
#  compute_flap_lift
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_interference_drag_highlift(wing,flapCDp,slatCDp):
    """Computes interference drag increment due to high-lift devices

    Assumptions:

    Source:
    Unknown

    Inputs:
    t_c                 (wing thickness ratio)                 [Unitless]
    flap_type                                                  [string]
    flap_c_chord        (flap chord as fraction of wing chord) [Unitless]
    flap_angle          (flap deflection)                      [radians]
    sweep               (Wing sweep angle)                     [radians]
    wing_Sref           (Wing reference area)                  [m^2]
    wing_affected_area  (Wing area affected by flaps)          [m^2]

    Outputs:
    dcl_max_flaps       (Lift coefficient increase)            [Unitless]

    Properties Used:
    N/A
    """          
    # Unpack inputs
    flap_type = wing.flaps.type 
    slat_type = wing.slats.type
    

    # Estimage kd factor
    if flap_type == 'plain:':
        Kintf = 0.0
    elif flap_type == 'single_slotted':
        Kintf = 0.40
    elif flap_type == 'single_slotted_Fowler' or flap_type == 'double_slotted_fixed_vane':
        Kintf = 0.25
    elif 'double_slotted' or flap_type == 'double_slotted_Fowler' or flap_type == 'triple_slotted_Fowler':
        Kintf = 0.25
    else:
        Kintf = 0.0

    if slat_type != 'None':
        Kints = 0.1
    else:
        Kints = 0

    dCDf_interf = Kintf * flapCDp   
    dCDs_interf = Kints * slatCDp

    return dCDf_interf, dCDs_interf 

