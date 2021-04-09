## @ingroup Methods-Aerodynamics-Fidelity_Zero-Drag
# compute_spoiler_lift.py
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

from SUAVE.Methods.Geometry.Three_Dimensional.compute_chord_length_from_span_location import compute_chord_length_from_span_location

# ----------------------------------------------------------------------
#  compute_spoiler_lift
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_spoiler_lift(vehicle,CLc):
    """Computes spoiler lift increment 

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
    bsi    = wing.spoilers.span_start
    bso    = wing.spoilers.span_end
    b      = wing.spans.projected
      
    bs  = (bso - bsi) / b
    
    dCL_sp = -CLc * bs / b
           
  
    return dCL_sp
