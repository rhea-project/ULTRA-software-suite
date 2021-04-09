## @ingroup Methods-Aerodynamics-Xfoil
#initialize_inputs.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

# SUAVE Imports
from SUAVE.Core import Data
# SUAVE-Xfoil Imports
from .create_xfoil_datastructure import create_xfoil_datastructure
from .write_geometry             import write_geometry

## @ingroup Methods-Aerodynamics-Xfoil
def initialize_inputs(geometry,M_eff,Cl_eff,Re_eff,ref_line,wing):
        """ This intializes the functions used in the AVL class

        Assumptions:
                None

        Source:
                None

        Inputs:
                xfoil_inputs - passed into the write_geometry, write_run_cases and write_input_deck functions

        Outputs:
                xfoil_inputs

        Properties Used:
                N/A
        """

        xfoil_inputs = create_xfoil_datastructure(geometry,ref_line,wing)
        write_geometry(xfoil_inputs)

        return xfoil_inputs
    
