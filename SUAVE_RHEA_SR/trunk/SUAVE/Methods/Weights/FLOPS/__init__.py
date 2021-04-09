## @defgroup Methods-Weights-FLOPS
# Provides structural weight correlations using the FLOPS method
# @ingroup Methods-Weights-FLOPS

""" SUAVE.Methods.Weights.FLOPS
    contains some useful methods or attributes
    for estimating weights with zero-order correlations
"""

# Attributes
from .                                import Operation_items
from .                                import Propulsion_systems
from .                                import System_and_equipment_items
from .empty                           import empty
from .systems                         import systems
from .weight_canard                   import weight_canard
from .weight_fin                      import weight_fin
from .weight_fuselage                 import weight_fuselage
from .weight_horizontal_tail          import weight_ht
from .weight_landing_gear   	      import weight_lg
from .weight_nacelles_air_induciton   import weight_nac
from .weight_paint                    import weight_paint
from .weight_payload_items            import weight_payload
from .weight_vertical_tail            import weight_vt
from .weight_wing                     import weight_wing

