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
def compute_parasite_drag_highlift(aero_clean,wing,alphad0,c_pr_c,Sfw,Ssw,Cla):
    """Computes induced drag increment due to high-lift devices

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

    # Unpack and prepare geometry
    Sref      = wing.areas.reference 
    bfi       = wing.flaps.span_start 
    bfo       = wing.flaps.span_end  
    cf        = wing.flaps.chord
    cs        = wing.slats.chord
    sweep     = wing.sweeps.quarter_chord
    flap_type = wing.flaps.type 
    df        = wing.flaps.angle  / Units.degrees

    if bfo - bfi != 0:
        if len(wing.Segments.keys())>0: 
            # obtain the geometry for each segment in a loop                                            
            n_segments = len(wing.Segments.keys())
            Cd0_sec    = np.zeros(n_segments)

            for i in range(n_segments):
                # Unpack airfoil data
                if wing.airfoils:
                    Cd0_sec[i] = wing.Segments[i].Airfoil.airfoil.Cd0 
                else: 
                    Cd0_sec[i] = 0.004

            Cd0 = np.average(Cd0_sec)

        else:
            # Unpack airfoil data
            Cd0 = wing.Airfoil.Cd0 

        # Estimage kd factor
        if flap_type == 'plain:':
            kd = 8.0178e-6*df**2 - 1.2614e-3*df + 2.6666e-1
        elif flap_type == 'single_slotted':
            kd = 1.7524e-7*df**4 - 2.8715e-5*df**3 + 1.6910e-3*df**2 - 4.0531e-2*df + 4.3149e-1
        elif flap_type == 'single_slotted_Fowler' or flap_type == 'double_slotted_fixed_vane':
            kd = 3.7766e-7*df**4 - 4.6409e-5*df**3 + 2.2632e-3*df**2 - 4.8505e-2*df + 4.7107e-1
        elif 'double_slotted' or flap_type == 'double_slotted_Fowler' or flap_type == 'triple_slotted_Fowler':
            kd = 7.2538e-8*df**4 - 1.4673e-5*df**3 + 1.0707e-3*df**2 - 3.2847e-2*df + 4.3584e-1
        else:
            kd = 0.0

        df    = df * Units.degrees
        dfCdp = kd * Cla * alphad0 * cf * df * np.sin(df) + Cd0 * (c_pr_c - 1)

        dfCDp = 1.15 * Sfw/(0.5*Sref) * dfCdp * np.cos(sweep)
       # dsCDp = aero_clean.wing[0,0,0,0] * Ssw/(0.5*Sref) * cs * np.cos(sweep)
        dsCDp = aero_clean.wing[0] * Ssw/(0.5*Sref) * cs * np.cos(sweep)

    else:
        dfCDp = 0.0
        dsCDp = 0.0
    
    return dfCDp, dsCDp
