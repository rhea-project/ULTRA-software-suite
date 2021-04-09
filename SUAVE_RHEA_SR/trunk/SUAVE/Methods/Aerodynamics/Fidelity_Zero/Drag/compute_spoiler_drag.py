## @ingroup Methods-Aerodynamics-Fidelity_Zero-Drag
# compute_spoiler_drag.py
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

from SUAVE.Methods.Geometry.Three_Dimensional import compute_chord_length_from_span_location

# ----------------------------------------------------------------------
#  compute_spoiler_drag
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_spoiler_drag(vehicle):
    """Computes spoiler drag increment 

    Assumptions:
        
    Source:
    Sadraey 'Aircraft Design'

    Inputs:


    Outputs:


    Properties Used:
    N/A
    """          

    # Unpack inputs
    Sref   = vehicle.reference_area

    wing   = vehicle.wings['main_wing']
    b      = wing.spans.projected
    Cs     = wing.spoilers.chord
    bsi    = wing.spoilers.span_start
    bso    = wing.spoilers.span_end
    delta  = wing.spoilers.angle   

    # Compute the spoiler area
    Csi = Cs * compute_chord_length_from_span_location(wing,0.5*bsi*b)        
    Cso = Cs * compute_chord_length_from_span_location(wing,0.5*bso*b)         
    bs  = (bso - bsi) * b
    Ssp = 0.5 * (Csi + Cso) * bs
    
    dCD_sp = 1.9 * np.sin(delta) * Ssp / Sref
           
  
    return dCD_sp
