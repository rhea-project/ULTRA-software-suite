## @defgroup Analyses-Aerodynamics Aerodynamics
# These are the analyses that control aerodynamic evaluations.
# @ingroup Analyses


from .Aerodynamics     import Aerodynamics
from .AVL              import AVL
from .AVL_Inviscid     import AVL_Inviscid 
from .Fidelity_Zero    import Fidelity_Zero
from .Markup           import Markup
from .Process_Geometry import Process_Geometry
from .Results          import Results
from .Supersonic_Zero  import Supersonic_Zero
from .Vortex_Lattice   import Vortex_Lattice
from .AERODAS          import AERODAS
from .Supersonic_OpenVSP_Wave_Drag import Supersonic_OpenVSP_Wave_Drag
from .Lifting_Line import Lifting_Line
from .Fidelity_Zero_Wing_Surrogate import Fidelity_Zero_Wing_Surrogate

#from .SU2_Euler        import SU2_Euler
#from .SU2_inviscid     import SU2_inviscid
#from .SU2_Euler_Super  import SU2_Euler_Super
#from .SU2_inviscid_Super import SU2_inviscid_Super

#from .SU2_Euler_wavedrag import SU2_Euler_wavedrag
#from .SU2_wave_drag       import SU2_wave_drag
