## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import pandas as pd
import math as m
from SUAVE.Core import Units
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#  Engine Weight Calculation
# ----------------------------------------------------------------------

def weight_engine(ac_type,NEW,NEF,THRUST,THRSO,WENGB,EEXP,inlet_included, \
                  nozzle_included,EINL,ENOZ,WINLB,WNOZB):
    
    """ Calculate the engine weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type 
        NEW             - number of wing mounted engines
        NEF             - number of fuselage mounted engines.                                                            The default value is 0 .
        THRUST          - Rated thrust of each scaled engine                                                      [N]
        THRSO           - Rated thrust of each baseline engine                                                    [N]
        WENGB           - Weight of baseline engine                                                               [kg]
        EEXP            - Engine weight scaling parameter. The default value is 1.15
        inlet_included  - equal 1 if the inlet weight is included in the baseline engine, otherwise 0.
        nozzle_included - equal 1 if the inlet weight is included in the baseline engine, otherwise 0.
        EINL            - Engine inlet weight scaling exponent. The default value is 1.
        ENOZ            - Engine nozzle weight scaling exponent. The default value is 1.
        WINLB           - Inlet weight for baseline engine                                                        [kg]
        WNOZB           - Nozzle weight for baseline engine                                                       [kg]
    Outputs:
        WENG = Weight of each scaled engine                                                                       [kg]                                             
      
    Properties Used:
        N/A
    """ 
    # unpack inputs
    THRUST = THRUST*0.22481                       #N to lb
    THRSO  = THRSO*0.22481                        #N to lb
    WENGB  = WENGB/Units.lb                       #kg to lb
    WINLB  = WINLB/Units.lb                       #kg to lb
    WNOZB  = WNOZB/Units.lb                       #kg to lb
    
    #Calculate the total number of engines
    NENG = NEW+NEF
    
    #Calculate the total number of nacelles
    TNAC = NENG+0.5*(NENG-2*m.floor(NENG/2))

    if WENGB==0:
        if ac_type=='transport' or ac_type=='HWB':
            WENGB=THRSO/5.5
        elif ac_type=='fighter' :
            WENGB=THRSO/8
        elif ac_type=='general aviation' :
            WENGB=THRSO/10.5
        else:
            print('Undefined aircraft type. Check the aircraft type definition')

    if inlet_included==0:
        WINL = WINLB*(THRUST/THRSO)**EINL
        
    if nozzle_included==0:
        WNOZL = WNOZB*(THRUST/THRSO)**ENOZ

    if EEXP>=0.3:  
        WENGP = WENGB*(THRUST/THRSO)**EEXP
    else:
        WENGP = WENGB+(THRUST-THRSO)*EEXP
        
    if (inlet_included==0 and nozzle_included==0) :   
        WENG  = WENGP + WINL + WNOZL
    else:
        WENG  = WENGP
    
    WENG = WENG * TNAC    
    WENG = WENG * Units.lb

    return WENG
      



