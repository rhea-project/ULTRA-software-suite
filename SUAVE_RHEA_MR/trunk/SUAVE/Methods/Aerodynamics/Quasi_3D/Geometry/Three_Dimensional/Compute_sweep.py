## @ingroup Methods-Geometry-Three_Dimensional
# AVL.py
#
# Created: Jul 2019, S. Karpuk 

# ----------------------------------------------------------------------
#  Imports
from math import atan as atan
from math import tan as tan

def Compute_sweep(b,Cr,Ct,c,sweep,tag):
    """The function computes a leading edge, a quarter-chord of a half-chord
       sweeps of the wing depending on the tag
       Transformations include:
           1. From leading edge to the specified sweep type
           2. From a specified sweep type to the leading edge 


    Assumptions:
    None

    Source:
    N/A

    Inputs:
    b                    [m]
    sweep                [rad]
    tag                  string

    Outputs:
    transformed_sweep    [rad]

    Properties Used:
    N/A
    """     
    
    if tag == 'from leading':
        transformed_sweep = atan((b*tan(sweep)-c*(Cr-Ct))/b)
    elif tag == 'to leading':
        transformed_sweep = atan((b*tan(sweep)+c*(Cr-Ct))/b)
    else:
        raise Exception('Check the input tag')
        

    return transformed_sweep
