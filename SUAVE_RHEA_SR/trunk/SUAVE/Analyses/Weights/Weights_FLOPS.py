## @ingroup Analyses-Weights
# Weights_FLOPS.py
#
# Created:  Dec 2019, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units,Data
from .Weights import Weights


# ----------------------------------------------------------------------
#  Analysis
# ----------------------------------------------------------------------

## @ingroup Analyses-Weights
class Weights_FLOPS(Weights):
    """ This is class that evaluates the weight of an aircraft using the FLOPS Method
    
    Assumptions:
        None

    Source:
        N/A

    Inputs:
        None
      
    Outputs:
        None
        
    Properties Used:
         N/A
    """
    def __defaults__(self):
        """This sets the default values and methods for the FLOPS
        aircraft weight analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """           
        self.tag = 'weights_FLOPS'
        
        self.vehicle  = Data()
        self.settings = Data()
        self.settings.weight_reduction_factors = Data()

        # Reduction factors are proportional (.1 is a 10% weight reduction)
        self.settings.weight_reduction_factors.main_wing = 0.                               # Main wing reduction factor
        self.settings.weight_reduction_factors.fuselage  = 0.                               # Fuselage reduction factor
        self.settings.weight_reduction_factors.empennage = 0. 		                    # Horizontal and vertical stabilizers reduction factor
        
        self.settings.load_factor = 1.5                                                     # Factor for GLA and MLA, 2.5 baseline and 1.0 1-g wing        
        
        self.settings.wing_weight_method = 1
        self.settings.wing_span_method   = 1
             
        
        
    def evaluate(self,conditions=None):
        """Evaluate the weight analysis.
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        results

        Properties Used:
        N/A
        """         
        # unpack
        vehicle  = self.vehicle
        settings = self.settings
        
        vehicle.settings.wing_weight_method = self.settings.wing_weight_method
        vehicle.settings.wing_span_method   = self.settings.wing_span_method

        empty    = SUAVE.Methods.Weights.FLOPS.empty

        
        # evaluate
        results = empty(vehicle)
        
        # storing weigth breakdown into vehicle
        vehicle.weight_breakdown = results 

        # updating empty weight
        vehicle.mass_properties.operating_empty = results.empty
              
        # done!
        return results

    '''def initialize(self):
        """Initializes the surrogate needed for lift calculation.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        self.geometry
        """                  
        self.process.compute.lift.inviscid_wings.geometry = self.geometry
        self.process.compute.lift.inviscid_wings.initialize()'''
           
    
    def finalize(self):
        """Finalize the weight analysis.
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """           
        self.mass_properties = self.vehicle.mass_properties
        
        return
        
