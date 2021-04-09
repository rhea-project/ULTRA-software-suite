## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
# aircraft_total.py
# 
# Created:  Dec 2013, A. Variyar,
# Modified: Feb 2014, A. Variyar, T. Lukaczyk, T. Orra 
#           Apr 2014, A. Variyar   
#           Jan 2016, E. Botero       


# ----------------------------------------------------------------------
#  Aircraft Total
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
def aircraft_total(state,settings,geometry):
    """Returns total aircraft lift and stores values

    Assumptions:
    None

    Source:
    None

    Inputs:
    state.conditions.aerodynamics.lift_coefficient    [Unitless]

    Outputs:
    aircraft_lift_total (lift coefficient)            [Unitless]

    Properties Used:
    N/A
    """      
    
    # unpack inputs
    conditions     = state.conditions

    spoiler_lift  = conditions.aerodynamics.spoiler_lift 

    aircraft_lift_total = state.conditions.aerodynamics.lift_coefficient + spoiler_lift

    return aircraft_lift_total
