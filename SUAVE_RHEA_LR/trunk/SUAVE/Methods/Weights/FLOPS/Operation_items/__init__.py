## @defgroup Methods-Weights-FLOPS-Operaton_items
# Provides structural weight correlations using the FLOPS method
# @ingroup Methods-Weights-FLOPS-Operation_items

""" SUAVE.Methods.Weights.FLOPS.Operation_items
    contains some useful methods or attributes
    for estimating weights with zero-order correlations
"""

# Attributes
from .cargo_containers        import weight_cargo_containers
from .crew_baggage            import weight_crew_baggage
from .engine_oil              import weight_engine_oil
from .fuel_capacity           import weight_fuel  
from .passenger_service       import weight_pass_service
from .unusable_fuel           import weight_unusable_fuel


