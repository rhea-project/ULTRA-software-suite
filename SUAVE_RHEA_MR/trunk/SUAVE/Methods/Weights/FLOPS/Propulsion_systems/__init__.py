## @defgroup Methods-Weights-FLOPS-Propulsion_systems
# Provides structural weight correlations using the FLOPS method
# @ingroup Methods-Weights-FLOPS-Propulsion_systems

""" SUAVE.Methods.Weights.FLOPS.Propulsion_systems
    contains some useful methods or attributes
    for estimating weights with zero-order correlations
"""

# Attributes
from .engine_weight                     import weight_engine
from .fuel_systems                      import weight_fuel_system
from .miscellaneous_propulsion_systems  import weight_engine_miscellaneous
from .thrust_reversers                  import weight_thrust_reversers



