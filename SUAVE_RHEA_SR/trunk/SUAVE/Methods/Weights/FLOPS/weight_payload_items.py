## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pandas as pd
from SUAVE.Core import Units
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#   Payload items Weight Calculation
# ----------------------------------------------------------------------

def weight_payload(ac_type,design_range,NPF,NPB,NPT,CARGOW,CARGOF,FULAUX):
    
    """ Calculate the weight of the payload items based on the FLOPS methods
       
    Source: 
        N/A
        
    Inputs:
      ac_type - 
      NPF     - Number of first class passengers
      NPB     - Number of business class passengers
      NPT     - Number of tourist class passengers
     CARGOW   - Cargo carried in wing  #including wing-mounted external stores for fighter/attack aircraft
     CARGOF   - Cargo carried in the fueselage  #including fuselage external sotres for fighter
                                              #cargo other than passenger baggage for trainsport aircraft
     FULAUX   - Auxiliary (external) fuel tank capacity for fighter/attack aircraft                       [N]

    Outputs:
      weight_payload - Total payload items weight                                                         [N]
                                                 
        
    Properties Used:
        N/A
    """
    
    # unpack inputs
    FULAUX       = FULAUX/Units.lb              # Convert kg to lb
    CARGOW       = CARGOW/Units.lb              # Convert kg to lb
    CARGOF       = CARGOF/Units.lb              # Convert kg to lb
    design_range = design_range/Units.nm        # Convert m to nm

    #Calculate the payload items weight using FLOPS methods      
    if ac_type=='fighter':
        WPBAG          = 1.15*FULAUX
        weight_payload = WPBAG+CARGOW+CARGOF

    else:
        NPASS  = NPF+NPB+NPT
        WPPASS = 165            #the default weight per passenger is 165lb
        WPASS  = WPPASS*NPASS
        if  1 <= design_range <= 900:
            BPP = 35
        elif 900 < design_range <= 2900:    
            BPP = 40
        else :
            BPP = 44
            
        WPBAG          = BPP*NPASS
        weight_payload = WPBAG+WPASS+CARGOW+CARGOF
        
    weight_payload = weight_payload*Units.lb
    
    return weight_payload

