## @ingroup Methods-Missions-Segments
# converge_root.py
# 
# Created:  Jul 2014, SUAVE Team
# Modified: Jan 2016, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import scipy.optimize
import numpy as np

from SUAVE.Core.Arrays import array_type

# ----------------------------------------------------------------------
#  Converge Root
# ----------------------------------------------------------------------

## @ingroup Methods-Missions-Segments
def converge_root(segment):
    """Interfaces the mission to a numerical solver. The solver may be changed by using root_finder.

    Assumptions:
    N/A

    Source:
    N/A

    Inputs:
    segment                            [Data]
    segment.settings.root_finder       [Data]
    state.numerics.tolerance_solution  [Unitless]

    Outputs:
    state.unknowns                     [Any]
    segment.state.numerics.converged   [Unitless]

    Properties Used:
    N/A
    """       
    
    unknowns = segment.state.unknowns.pack_array()
    
    # Find the normalization factors
    segment.state.unknowns_normalization_factor  = unknowns*1.
    segment.state.unknowns_normalization_factor[segment.state.unknowns_normalization_factor==0] = 1e-16 
    
    # Run one iteration to get the scaling
    segment.process.iterate(segment)
    segment.state.residual_normalization_factor = 1*segment.state.residuals.pack_array() 
    segment.state.residual_normalization_factor[segment.state.residual_normalization_factor==0] = 1e-16
    
    # Normalize the unknowns
    unknowns = unknowns/segment.state.unknowns_normalization_factor
    
    
    try:
        root_finder = segment.settings.root_finder
    except AttributeError:
        root_finder = scipy.optimize.fsolve 
     
    unknowns,infodict,ier,msg = root_finder( iterate,
                                          unknowns,
                                          args = segment,
                                          xtol = segment.state.numerics.tolerance_solution,
                                          full_output=1)
    if ier!=1:
        print("Segment did not converge. Segment Tag: " + segment.tag)
        print("Error Message:\n" + msg)
        segment.state.numerics.converged = False
    else:
        segment.state.numerics.converged = True
    
    '''ier = 0   
    if ier!=1:
        print("solving segment with lm")
        OptimizeResult = scipy.optimize.root( iterate,unknowns,
                                        method = 'lm',
                                        args = segment,
                                        tol = segment.state.numerics.tolerance_solution)
        if not ier :
            segment.state.numerics.converged = False
            print(" segment did not converged with the lm solver: " + segment.tag )
        else:
            print(" segment converged with the lm solver " + segment.tag)
            segment.state.numerics.converged = True '''
         
                            
    return
    
# ----------------------------------------------------------------------
#  Helper Functions
# ----------------------------------------------------------------------

## @ingroup Methods-Missions-Segments
def iterate(unknowns, segment):
    
    """Runs one iteration of of all analyses for the mission.

    Assumptions:
    N/A

    Source:
    N/A

    Inputs:
    state.unknowns                [Data]
    segment.process.iterate       [Data]

    Outputs:
    residuals                     [Unitless]

    Properties Used:
    N/A
    """       
    if isinstance(unknowns,array_type):
        segment.state.unknowns.unpack_array(unknowns)
    else:
        segment.state.unknowns = unknowns
        
    segment.process.iterate(segment)
    
    residuals = segment.state.residuals.pack_array()
        
    return residuals 
