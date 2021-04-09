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
#   Paint Weight Calculation
# ----------------------------------------------------------------------

def weight_paint(area_dens_paint,wet_ara_wing,wet_ara_ht,wet_ara_vt, \
                 wet_ara_fues,wet_ara_nac,wet_ara_can,wet_ara_fin,NFUSE,NVTAIL):
    
    """ Calculate paint weight of the aircraft based on the FLOPS methods
            
    Inputs:
        area_dens_paint : Area density of paint for all wetted areas [kg/m²]                             
        wet_ara_wing    : Wetted area of wings                       [m²]
        wet_ara_ht      : Wetted area of horizontal tail             [m²]
        wet_ara_vt      : Wetted area of vertical tail               [m²]
        wet_ara_fues    : Wetted area of fueselage                   [m²]
        wet_ara_nac     : Wetted area of nacelles                    [m²]
        wet_ara_can     : Wetted area of canards                     [m²]
        wet_ara_fin     : Wetted area of fin                         [m²]
        NFUSE           : Number of fuselages
        NVTAIL          : Number of Vtail
    Outputs:
       WVT                                                           [kg]
        
    Properties Used:
        N/A
    """ 
     # unpack inputs
    area_dens_paint = area_dens_paint/(Units.lb/Units.ft**2)    # Convert from kg per meter squared to lb per ft squared 
    wet_ara_wing    = wet_ara_wing/Units.ft**2                  # Convert meters squared to ft squared
    wet_ara_ht      = wet_ara_ht/Units.ft**2                    # Convert meters squared to ft squared
    wet_ara_vt      = wet_ara_vt/Units.ft**2                    # Convert meters squared to ft squared
    wet_ara_fues    = wet_ara_fues/Units.ft**2                  # Convert meters squared to ft squared
    wet_ara_nac     = wet_ara_nac/Units.ft**2                   # Convert meters squared to ft squared
    wet_ara_can     = wet_ara_can/Units.ft**2                   # Convert meters squared to ft squared
    wet_ara_fin     = wet_ara_fin/Units.ft**2                   # Convert meters squared to ft squared

    wet_ara_fues    = wet_ara_fues * NFUSE
    wet_ara_vt      = wet_ara_vt * NVTAIL

    #Calculate weight of the paint using  FLOPS methods
    paint_weight = area_dens_paint*(wet_ara_wing+wet_ara_ht+wet_ara_vt+wet_ara_fues \
                                   +wet_ara_nac+wet_ara_can+wet_ara_fin)
    paint_weight = paint_weight*Units.lb

    
    return paint_weight
