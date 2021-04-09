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
#  Engine Weight Calculation
# ----------------------------------------------------------------------

def weight_fuel_system(NEW,NEF,ac_type,FMXTOT,VMAX,NTANK):
    
    """ Calculate the fuel system weight (including the fuel tanks and plumbing) weight  based on the FLOPS methods    
   
    Source: 
        N/A
        
    Inputs:
    
        NEW     - number of wing mounted engines
        NEF     - number of fuselage mounted engines. The default value is 0 .
        FMXTOT  - Aircraft maximum fuel capacity                                            [kg]
        ac_type - Aircraft type
        VMAX    - Maximum Mach number
        NTANK   - Number of tanks, for fighter aircraft
        
    Outputs:
        WPMSC = Total fuel system weight weight                                             [kg]                                             
      
    Properties Used:
        N/A
    """ 
    # unpack inputs
    FMXTOT = FMXTOT/Units.lb
    NENG   = NEW+NEF
    
    #Calculate the engine controls weight
    if ac_type=='transport':
        WFSYS = 1.07*FMXTOT**0.58*NENG**0.43*VMAX**0.34
    elif ac_type=='fighter':
        WFSYS = 36*FMXTOT**0.2*NTANK**0.5*NENG**0.4
    elif ac_type=='general aviation':
        WFSYS = 1.07*FMXTOT**0.58*NENG**0.43
        
    WFSYS=WFSYS*Units.lb
    
    return WFSYS
      


