## @defgroup Methods-Aerodynamics-AVL_Loads AVL_Loads
# Functions to AVL Loads calculations
# @ingroup Methods-Aerodynamics

""" SUAVE AVL Interface Package Setup
"""

from .write_input_deck         import write_input_deck
from .write_run_cases          import write_run_cases
from .run_analysis             import run_analysis
from .write_geometry           import write_geometry
from .create_avl_datastructure import create_avl_datastructure
from .translate_data           import translate_conditions_to_cases, translate_results_to_conditions

from . import Data
