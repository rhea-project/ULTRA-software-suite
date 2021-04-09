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
def compute_induced_drag_highlift(wing,dfCL0):
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
    AR    = wing.aspect_ratio   
    bfi   = wing.flaps.span_start 
    bfo   = wing.flaps.span_end  
    sweep = wing.sweeps.quarter_chord

    # interpolate data
    bfi_data = np.array([0.0, 0.2, 0.4, 0.6])
       
    for i in range(len(bfi_data)):
        if bfi >= bfi_data[i]:
            bf_interp = bfi_data[i]

    Kv_data = np.zeros(2)
    if bf_interp == 0.0:
        b_range    = np.array([bf_interp, bf_interp + 0.2])
        a1         = np.array([-0.4133, 0.1936, -0.9886, 0.463, -0.3324, -0.3137, 0.3906])
        Kv_data[0] = a1[0]*bfo**a1[1] + a1[2]*bfo**a1[3]*AR**a1[4] + AR**a1[5] + a1[6]
        a2         = np.array([-3.095, 0.03409, -2.022, 0.04482, 0.08923, 0.1476, 4.301])
        Kv_data[1] = a2[0]*bfo**a2[1] + a2[2]*bfo**a2[3]*AR**a2[4] + AR**a2[5] + a2[6]
    elif bf_interp == 0.2:
        b_range    = np.array([bf_interp, bf_interp + 0.2])
        a1         = np.array([-3.095, 0.03409, -2.022, 0.04482, 0.08923, 0.1476, 4.301])
        Kv_data[0] = a1[0]*bfo**a1[1] + a1[2]*bfo**a1[3]*AR**a1[4] + AR**a1[5] + a1[6]
        a2         = np.array([3.119e-05, -14.52, -0.1112, 1.229, 0.2547,  0.446, -0.02496])
        Kv_data[1] = a2[0]*bfo**a2[1] + a2[2]*bfo**a2[3]*AR**a2[4] + a2[5]*AR**a2[6] 
    elif bf_interp == 0.4:
        b_range    = np.array([bf_interp, bf_interp + 0.2])
        a1         = np.array([3.119e-05, -14.52, -0.1112, 1.229, 0.2547,  0.446, -0.02496])
        Kv_data[0] = a1[0]*bfo**a1[1] + a1[2]*bfo**a1[3]*AR**a2[4] + a1[5]*AR**a1[6] 
        a2         = np.array([4.081, -3.55, 1.219, 0.02192])
        Kv_data[1] = a2[0]*bfo + a2[1]*bfo**a2[2]*AR**a2[3]   
    elif bf_interp == 0.6:
        b_range    = np.array([bf_interp, bf_interp])
        a1         = np.array([4.081, -3.55, 1.219, 0.02192])
        Kv_data[0] = a1[0]*bfo + a1[1]*bfo**a1[2]*AR**a1[3] 
        Kv_data[1] = a1[0]*bfo + a1[1]*bfo**a1[2]*AR**a1[3]   
    else:
        b_range    = np.array([bf_interp, bf_interp])
        Kv_data[0] = 0.0;    Kv_data[1] = 0.0

    Kinv = np.interp(bfi,b_range,Kv_data)

    dfCdi = Kinv**2 * dfCL0**2 * np.cos(sweep)
    

    return dfCdi
