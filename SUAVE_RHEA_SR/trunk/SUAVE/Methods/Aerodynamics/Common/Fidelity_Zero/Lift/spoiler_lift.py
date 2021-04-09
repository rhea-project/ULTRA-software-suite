## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
# spoiler_lift.py
#
# Created:  

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
#  Adds the spoiler lift during the main mission
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Common-Fidelity_Zero-Lift
def spoiler_lift(state,settings,geometry):
    """Adds a spoiler drag increment

    Assumptions:
    None

    Source:
    None

    Inputs:
    settings.spoiler_drag_increment  [Unitless]

    Outputs:
    spoiler_drag                     [Unitless]

    Properties Used:
    N/A
    """    
    
    # unpack inputs
    conditions     = state.conditions
    configuration  = settings
    drag_breakdown = conditions.aerodynamics.drag_breakdown

    # various drag components
    spoiler_lift = settings.spoiler_lift_increment

    # untrimmed drag
    conditions.aerodynamics.spoiler_lift = spoiler_lift
    
    return spoiler_lift
