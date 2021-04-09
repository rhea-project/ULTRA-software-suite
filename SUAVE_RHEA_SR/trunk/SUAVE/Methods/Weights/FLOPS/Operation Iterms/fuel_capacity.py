## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------


import pandas as pd
from SUAVE.Core import Units

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Total fuel Weight Calculation
# ----------------------------------------------------------------------

def weight_fuel(ac_type,
                FULDEN,
                FWMAX,
                SW,
                TCA,
                TR,
                SPAN,
                FULWMX,
                FUELM,
                FUFU):
    """ Calculate the fuel weight  based on the FLOPS methods
    
   
    Source: 
        N/A
        
    Inputs:
        ac_type - aircraft type                 : 1 for BWB
                                                  2 for fighter/attack aircraft
        FULDEN   - Fuel density ratio for alternate fuels compared to jet fuel  #default value is 1
        FWMAX    - Factor for wing fuel capacity             #default value is 23
        SW       - Reference wing area                       [m²]
        TCA      - Weighted average of the wing thickness to chord ratio
        TR       - Taper ratio of the wing
        SPAN     - Wing span                                 [m]
        FULWMX   - Total fuel capacity of the wing           [kg]
        FUELM    - Total aircraft fuel weight                [kg]
        FUFU     - Total fuel capacity of the fuselage       [kg]
        
    Outputs:
      FMXTOT     - Aircraft maximum fuel capacity            [kg]                                             
      
    Properties Used:
        N/A
    """ 
     # unpack inputs
    SW=SW/Units.ft**2 # Convert meters squared to ft squared
    SPAN=SPAN/Units.ft
    FULWMX=FULWMX/Units.lb # Convert kg to lb
    FUELM=FUELM/Units.lb
    FUFU=FUFU/Units.lb
    #Wing fuel capacity
    if FULWMX==0:
       FULWMX= FULDEN*FWMAX*SW**2*TCA*(1-(TR/((1+TR)**2)))/SPAN
       if FUELM==0:
           FUELM=FULWMX
    if FUELM!=0 :
       if  FUELM>FULWMX and FUFU==0:
           FUFU=FUELM-FULWMX
       elif FUFU!=0 and FUELM>FULWMX+FUFU :
            FULAUX=FUELM-(FULWMX+FUFU) if ac_type==1 else 0
   
    FMXTOT=FULWMX+FUFU
    FMXTOT=FMXTOT*Units.lb
    FULWMX=FULWMX*Units.lb
    FUFU=FUFU*Units.lb # Convert lb to kg
    print(FULWMX,FUFU)
    return FMXTOT
      
 