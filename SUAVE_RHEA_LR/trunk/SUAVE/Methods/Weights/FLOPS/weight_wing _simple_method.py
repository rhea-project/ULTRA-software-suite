## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import math as m
from SUAVE.Core import Units
import pandas as pd
import matplotlib.pyplot as plt
# ----------------------------------------------------------------------
#   Wing Weight Calculation
# ----------------------------------------------------------------------

def weight_wing(DG,
                WSR,
                S_wing,
                wing_sweep,
                wing_TR,
                wing_TCA,
                wing_span,
                wing_AR,
                wing_span_method=1,
                aircraft_type="transport",
                WF=5,
                OSSPAN=20,
                FSTRT=0,
                FAERT=1):
    """ Calculate the wing weight of the aircraft based on the FLOPS methods
    
    Source: 
        N/A
        
    Inputs:
        DG   -Design gross weight                     [kg] 
        WSR - Required wing loading                   [kg/m**2]
        S_wing - Reference wing area                  [m**2]
        wing_sweep- sweep of the wing                 [deg]
        wing_TR - wing Taper ratio                    [] 
        wing_TCA - Weighted average of the wing thickness to chord ratio # if not available default value 12%
        wing_span - span of the wing                  [m]
        wing_AR - Wing Aspect Ratio 
        GLOV - Total glove and bat area beyond theoretical wing area [m**2]
        wing_span_method                              # default value 1 for  fixed wing span method
        ETA - Local wing station location             [m]  # used for detailed wing method, default value 0
        chord_length   - Local chord length           [m]  # used for detailed wing method, default value 0
        aircraft_type                                 # default 'transport', otherwise 'HWB' for hybrid wing
        WF - Maximum fuselage width                   [m] #used only for HWB aircraft, deufal value 5
        OSSPAN - Outboard wing semispan of HWB aircraft [m] #used only for HWB aircraft, deufal value 20
        FAERT - Aeroelastic tailoring factor used in the design of the wing : 0.0 for no aeroelastic tailoring 
                                                                              1.0 for maximum aeroelastic tailoring 
        FSTRT: Wing strut bracing factor                                    : 0.0 for no wing strut  
                                                                              1.0 for full benefit from strut bracing
        
        
       
    Outputs:
      weight_wing         Weight of the wing using the FLOPS simple method                        
        
    Properties Used:
        N/A
    """ 
    # unpack inputs
    DG=DG/Units.lb # Convert kg to lb
    WSR=WSR/Units.lb/Units.ft**2 # Convert from kg per meter squared to lb per ft squared
    S_wing=S_wing/Units.ft**2
    wing_span=wing_span/Units.ft #convert meter to ft 
    wing_sweep=m.radians(wing_sweep) #convert grad to radiant
  

    #Calculate weight of wing using FLOPS methods   
      #Simplified Wing Weight Estimation Method
    wing_EMS=1        #Calculate Wing strut bracing factor
    if wing_AR<=5: 
        wing_CAYA=0
    else: 
        wing_CAYA= wing_AR-5
    wing_TLAM=m.tan(wing_sweep)-(2*(1-wing_TR))/(wing_AR*(1+wing_TR))
    wing_SLAM=wing_TLAM/m.sqrt(1+wing_TLAM)
    C_4=1-0.5*FAERT
    C_6=0.5*FAERT-0.16*FSTRT
    wing_CAYL=(1-wing_SLAM**2)*(1+C_6*wing_SLAM**2+0.03*wing_CAYA*C_4*wing_SLAM) #Calculate  The Wing sweep factor
    wing_eq_mat_fac=(0.215*(0.37+0.7*wing_TR)*m.pow((wing_span**2/S_wing),wing_EMS))/(wing_CAYL*wing_TCA)  #Calculate The equivalent bending factor
    W1N1R=8.8*wing_eq_mat_fac*(1+m.sqrt(6.25/wing_span))*3.75*wing_span*(1-0.1*FAERT)*(0.5/10**6) # wing bending material weight without the efects of inertia relief
    W2=0.68*(0.25*S_wing)**(0.34)*DG**(0.6)      #movable wing surface is assumed to be 25% of the total surface                                                                                          #Total Wing Shear Material and Control Surface Weight
    W3=0.035*S_wing**(1.5)
    CAYE=1-0.03*1  # 1 for single wing mounted engine
    W1=((DG*CAYE*W1N1R+W2+W3)/(1+W1N1R))-W2-W3
    print(W1,W2,W3,'-----')
    weight_wing=W1+W2+W3
    weight_wing=weight_wing*Units.lb
    return weight_wing
    

