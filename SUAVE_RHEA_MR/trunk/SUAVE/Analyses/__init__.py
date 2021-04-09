## @defgroup Analyses
# Analyses are classes that are comprised of several methods. 
# Methods operate on attributes. This process and structure is described 
# <a href="http://suave.stanford.edu">here</a>.

from .Analysis  import Analysis
from .Sizing    import Sizing
from .Surrogate import Surrogate
from .Process   import Process
from .Settings  import Settings
from .Vehicle   import Vehicle

from . import Aerodynamics
from . import Aeroelasticity
from . import Stability
from . import Energy
from . import Weights
from . import Geometry
from . import Loads
from . import Mission
from . import Structures
from . import Atmospheric
from . import Planets
from . import Sizing
from . import Noise
from . import Costs
from . import Constraint_Analysis
