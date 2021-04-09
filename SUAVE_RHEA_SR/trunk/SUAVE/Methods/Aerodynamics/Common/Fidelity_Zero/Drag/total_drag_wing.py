## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Drag
# parasite_drag_wing.py
# 
# Created:  Dec 2013, SUAVE Team
# Modified: Jan 2016, E. Botero       

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# local imports
from SUAVE.Methods.Aerodynamics.Common.Fidelity_Zero.Drag import parasite_drag_wing, induced_drag_aircraft, compressibility_drag_wing
from smt.surrogate_models import RMTB

# suave imports
from SUAVE.Core import Data, Units

# package imports
import numpy as np

# ----------------------------------------------------------------------
#   Parasite Drag Wing
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Drag
def total_drag_wing_surrogate(state,settings,geometry):
    """Computes the the wing drag based on available surrogate data

    Assumptions:
    

    Source:


    Inputs:
    settings.wing_parasite_drag_form_factor      [Unitless]
    state.conditions.freestream.
      mach_number                                [Unitless]
      temperature                                [K]
      reynolds_number                            [Unitless]
    geometry.
      areas.reference                            [m^2]
      chords.mean_aerodynamic                    [m]
      thickness_to_chord                         [Unitless]
      sweeps.quarter_chord                       [radians]
      aspect_ratio                               [Unitless]
      spans.projected                            [m]
      areas.exposed                              [m^2]
      areas.affected                             [m^2]
      areas.wetted                               [m^2]
      transition_x_upper                         [Unitless]
      transition_x_lower                         [Unitless]
      
      
    Outputs:
    wing_parasite_drag                           [Unitless]

    Properties Used:
    N/A
    """
    
    # unpack
    conditions             = state.conditions
    wings                  = geometry.wings
    fuselages              = geometry.fuselages
    propulsors             = geometry.propulsors
    vehicle_reference_area = geometry.reference_area

    CL   = conditions.aerodynamics.lift_coefficient
    CDp0 = conditions.aerodynamics.drag_breakdown.parasite.main_wing.parasite_drag_coefficient
    H   = conditions.freestream.altitude
    M   = conditions.freestream.mach_number

    CD_surrogate = settings.surrogate_drag_wing

    CRUD = wings['main_wing'].CRUD
    misc = wings['main_wing'].miscellaneous
    IF   = wings['main_wing'].interference 

    eps   = 1e-6
    error = 1

    while error > eps:
        # compute parasite drag 
        total_wing_drag = 0.0

        # main wing parasite drag from the surrogate model
        xy      = np.hstack([CL.reshape(-1, 1),H.reshape(-1, 1),M.reshape(-1, 1)])
        CD_wing = CD_surrogate.predict_values(xy)
        #conditions.aerodynamics.drag_breakdown.parasite['main_wing'].parasite_drag_coefficient = CD_wing * wings['main_wing'].areas.reference/vehicle_reference_area
        total_wing_drag += CD_wing * wings['main_wing'].areas.reference/vehicle_reference_area

        # dump to condtitions
        state.conditions.aerodynamics.drag_breakdown.parasite.total = total_wing_drag 

        # Compute compressibility drag
        CDc = compressibility_drag_wing(state,settings,wings['main_wing'])
        state.conditions.aerodynamics.drag_breakdown.compressible.total = CDc

        # Compute induced drag
        #CDi = induced_drag_aircraft(state,settings,geometry)
        CDi = CL**2/(np.pi*0.9*wings['main_wing'].aspect_ratio)

        CDp = total_wing_drag - CDi

        error = np.max(np.abs(CDp-CDp0))
        CDp0  = CDp

    total_parasite_drag = 0

    # Correct the final wing paraside drag with additional drag components
    total_parasite_drag += (CDp+misc) * IF * ( 1 + CRUD )

    # parasite drag from remaining components
    for wing in wings.values():
        if wing.tag != 'main_wing':
            parasite_drag = conditions.aerodynamics.drag_breakdown.parasite[wing.tag].parasite_drag_coefficient 
            conditions.aerodynamics.drag_breakdown.parasite[wing.tag].parasite_drag_coefficient = parasite_drag * wing.areas.reference/vehicle_reference_area
            total_parasite_drag += parasite_drag * wing.areas.reference/vehicle_reference_area

    # Finish parasite drag calculations for remaining aircraft components
    # from fuselage
    for fuselage in fuselages.values():
        if fuselage.tag == 'fuselage_bwb':
            continue
        parasite_drag = conditions.aerodynamics.drag_breakdown.parasite[fuselage.tag].parasite_drag_coefficient 
        conditions.aerodynamics.drag_breakdown.parasite[fuselage.tag].parasite_drag_coefficient = parasite_drag * fuselage.areas.front_projected/vehicle_reference_area
        total_parasite_drag += parasite_drag * fuselage.areas.front_projected/vehicle_reference_area
    
    # from propulsors
    for propulsor in propulsors.values():
        ref_area = propulsor.nacelle_diameter**2 / 4 * np.pi
        parasite_drag = conditions.aerodynamics.drag_breakdown.parasite[propulsor.tag].parasite_drag_coefficient 
        conditions.aerodynamics.drag_breakdown.parasite[propulsor.tag].parasite_drag_coefficient  = parasite_drag * ref_area/vehicle_reference_area * propulsor.number_of_engines
        total_parasite_drag += parasite_drag * ref_area/vehicle_reference_area * propulsor.number_of_engines
 
    # from pylons
    try:
        parasite_drag = conditions.aerodynamics.drag_breakdown.parasite['pylon'].parasite_drag_coefficient
    except:
        parasite_drag = 0. # not currently available for supersonics

    total_parasite_drag += parasite_drag
        
    # Additional friction drag penalties
    total_parasite_drag += geometry.misc  

    # dump to condtitions
    state.conditions.aerodynamics.drag_breakdown.parasite.total = total_parasite_drag 


    return total_parasite_drag
