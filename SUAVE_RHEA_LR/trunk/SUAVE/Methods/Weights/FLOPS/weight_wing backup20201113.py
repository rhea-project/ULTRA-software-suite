## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk


# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import math as m
from SUAVE.Core import Units
# ----------------------------------------------------------------------
#   Wing Weight Calculation
# ----------------------------------------------------------------------

def weight_wing(DGW,GW,WSR,S_wing,wing_sweep,wing_TR, wing_TCA,wing_span,wing_AR,GLOV,
                wing_span_method,wing_weight_method,aircraft_type,WF,OSSPAN,FSTRT,
                FAERT,PCTL,SFLAP,eng_num,FCOMP,FNEF,SAFTB,RSPCHD,RSPSOB,XLP,Nult,CAYF,
                cab_LE,VARSWP,ETA=[0],chord_length=[0]):
    
    """ Calculate the wing weight of the aircraft based on the FLOPS methods
    
    Source: 
        N/A
        
    Inputs:
        DGW            - Design gross weight                                                      [kg]
        GW             - Ramp weight                                                              [kg]
        WSR            - Required wing loading                                                    [kg/m**2]
        S_wing         - Reference wing area                                                      [m**2]
        wing_sweep     - sweep of the wing                                                        [radians]
        wing_TR        - wing Taper ratio                                                         [Unitless] 
        wing_TCA       - Weighted average of the wing thickness to chord ratio 
        wing_span      - span of the wing                                                         [m]
        wing_AR        - Wing Aspect Ratio 
        GLOV           - Total glove and bat area beyond theoretical wing area                    [m**2]
        wing_span_method                              1 for  fixed wing span method
        wing_weight_method                            1 for  simple method
                                                      2 for  detailed method
        aircraft_type                                                                                    # default 'transport', otherwise 'HWB' for hybrid wing
        WF             - Maximum fuselage width                                                   [m]    # used only for HWB aircraft, deufal value 5
        OSSPAN         - Outboard wing semispan of HWB aircraft                                   [m]    # used only for HWB aircraft, deufal value 20
        FSTRT          - Wing strut bracing factor : 0.0 for no wing strut  
                                                     1.0 for full benefit from strut bracing
        FAERT          - Aeroelastic tailoring factor used in the design of the wing  :
                         0.0 for no aeroelastic tailoring 
                         1.0 for maximum aeroelastic tailoring
        FCOMP          - Composite utilization
                         0.0 for no composites
                         1.0 for maximum use of composites
        PCTL           - Fraction of load carried by the defined wing
        ETA            - Local wing station location                                               [m]  # used for detailed wing method, default value 0
        SFLAP          - flap total area                                                           [m**2]
        FNEF           - number of fuselage-mounted engines  
        SAFTB          - aftbody area (used for the HWB)                                           [m**2]
        RSPCHD         - Percent chord of the HWB fuselage rear spar at the fuselage centerline.
                         The default value is 70%
        RSPSOB         - Percent chord of the HWB fuselage rear spar at the side of body.
                         By defaultRSP SOB=RSPCHD.  Where RSPCHD is the percent chord of the
                         HWB fuselage rear spar at the fuselage centerline
        XLP            - Passenger compartment length                                              [m]
        cab_LE         - passenger cabin leading-edge sweep                                        [radians]
        CAYF           - multiple fuselage factor (1 - single, 0.5 - multiple)
        Nult           - Ultimate load factor
        VARSWP         - Variable sweep factor (0 - fixed sweep, 1 - variable sweep)
        chord_length   - Local chord length                                                        [m]  # used for detailed wing method, default value 0
        
        
       
    Outputs:
      weight_wing                                                                                  [kg]                    
        
    Properties Used:
        N/A
        
    """
    
    # unpack inputs                  
    GW           = GW/Units.lb                    # Convert kg to lb
    WSR          = WSR/(Units.lb/Units.ft**2)     # Convert meters squared to ft squared 
    S_wing       = S_wing/Units.ft**2             # Convert from meter squared to ft squared
    wing_span    = wing_span/Units.ft             # Convert meter to ft
    WF           = WF/Units.ft
    OSSPAN       = OSSPAN/Units.ft
    SFLAP        = SFLAP/Units.ft**2
    GLOV         = GLOV/Units.ft**2
    SAFTB        = SAFTB/Units.ft**2
    XLP          = XLP/Units.ft 

    #Calculate weight of wing using FLOPS methods
    #Calculate Design Gross Weight
    if DGW>0 and  DGW<5:
       DG = DGW*GW
    elif DGW>5:
         DG = GW
    else:
         DG = GW
    #Calculate wing reference area
    if WSR!=0:
          S_wing = GW/WSR
    elif wing_weight_method==1 :
          S_wing = (wing_span**2)/wing_AR+GLOV
    elif wing_weight_method ==2:
         S_wing = 0
         for i in range(len(ETA)-1):
            S_wing += (ETA[i+1]-ETA[i])*(C[i]+C[i+1])
    #Calculate wing trapezoidal Area
    wing_trapez_area = S_wing- GLOV
    #Calculate Wing Span and Wing Aspect Ratio 
    if wing_span_method==1:
       wing_span_B = fixed_span = wing_span                                     # Calculate wing span
       wing_AR     = (fixed_span**2)/wing_trapez_area                           # Calculate wing aspect ratio
    else :
        if wing_AR!=0:
           if wing_span==0:
              if aircraft_type=="HWB":
                 wing_span = WF+OSSPAN*2
              else: 
                 wing_span = m.sqrt(wing_AR*(S_wing-GLOV))
           else:
              wing_AR = (wing_span_B**2)/(S_wing-GLOV)
              
        wing_span_B = m.sqrt(wing_AR*wing_trapez_area)
 
    if aircraft_type == 'transport' or aircraft_type == 'HWB':
        A = [8.8,6.25,0.68,0.34,0.6,0.035,1.5]
    elif aircraft_type == 'fighter':
        A = [6.8,0.0,0.12,0.65,0.62,0.8,1.2]
    else:
        A = [30.0,0.0,0.25,0.5,0.5,0.16,1.2]

        
    ''' 2020-11-02, for Strut-braced wing configuration, as follows
    #Simplified Wing Weight Estimation Method
    wing_EMS = 1-0.25*FSTRT                                                                                     # Calculate Wing strut bracing factor
    if wing_AR<=5: 
        wing_CAYA = 0
    else: 
        wing_CAYA = wing_AR-5
 
    wing_TLAM       = m.tan(wing_sweep)-(2*(1-wing_TR))/(wing_AR*(1+wing_TR))
    wing_SLAM       = wing_TLAM/m.sqrt(1+wing_TLAM**2)
    C_4             = 1-0.5*FAERT
    C_6             = 0.5*FAERT-0.16*FSTRT
    wing_CAYL       = (1-wing_SLAM**2)*(1+C_6*wing_SLAM**2+0.03*wing_CAYA*C_4*wing_SLAM)                    # Calculate  The Wing sweep factor
    wing_eq_mat_fac = 0.215*(0.37+0.7*wing_TR)*m.pow(wing_span**2/S_wing,wing_EMS)/(wing_CAYL*wing_TCA)    # Calculate The equivalent bending factor

    if VARSWP <= 0:
        VFACT = 1
    else:
        VFACT = 1+VARSWP*(0.96/m.cos(wing_sweep)-1)
        
    W1N1R           = A[0]*wing_eq_mat_fac*(1+m.sqrt(A[1]/wing_span))*Nult*wing_span*(1-0.4*FCOMP)* \
                      (1-0.1*FAERT)*CAYF*VFACT*(PCTL/1000000)                                                  # wing bending material weight without the effects of inertia relief
    W2              = A[2]*(1-0.17*FCOMP)*SFLAP**A[3]*DG**A[4]                                                             # Total Wing Shear Material and Control Surface Weight
    W3              = A[5]*(1-0.3*FCOMP)*S_wing**A[6]
    CAYE            = 1-0.03*eng_num                                                                        # eng_num for single wing mounted engine
    W1              = (DG*CAYE*W1N1R+W2+W3)/(1+W1N1R)-W2-W3

    if aircraft_type == 'HWB':
        XLW    = XLP-m.tan(cab_LE)*0.5*WF
        TRAFTB = ((1-RSPSOB)*XLW/RSPSOB)/((1-RSPCHD)*XL)
        W4     = (1+0.05*FNEF)*0.53*SAFTB*DG**0.2*(0.5+TRAFTB)*(1-0.17*FCOMP)
    else:
        W4 = 0
   
    weight_wing     = W1+W2+W3+W4 #+ 0.5*S_wing   # INCLUDES DOUBLE-SLOTTED FLAP RATIO
    weight_wing     = weight_wing*Units.lb'''




    # Wing and Strut weight estimation method for Strut-braced wing (SBW) aircraft
    '''# Source:  P. Chiozzotto G. Initial weight estimate of advanced transport aircraft concepts considering aeroelastic effects[C]//55th AIAA aerospace sciences meeting. 2017: 0009.

    # Parameters for exponents E and constant C (related to the aircraft configuration and materials), follows parameters valuse are for the SBW configuration with carbon fiber reforced plastic (CFRP)
    
    # Constant C
    C_c     =   0.00225     # Constant C for wing covers
    C_wr    =   0.209       # Constant C for wing ribs and webs
    C_sj    =   0.00101     # Constant C for wing strut and juries
    C_a     =   464         # for aileron efficiency

    # engine relief factor
    ke_c    =   0.99        # covers
    ke_wr   =   0.984       # webs and ribs
    ke_sj   =   0.945       # strut and juries
    
    # Exponents E for wing covers
    E_c_m     = 1.351       # m_TO
    E_c_w     = -0.708      # W/S
    E_c_A     = 1.19        # aspect ratio A
    E_c_sweep = -1.794      # wing sweep Lamda
    E_c_t     = -0.724      # t/c
    E_c_V     = 0.02        # flight speed V
    E_c_l     = 0.603       # taper ratio lamda
    E_c_nz    = 0.886       # nz design positive limit load factor
    E_c_e     = 1.511       # strut position eta
    
    # Exponents E for wing webs and ribs
    E_wr_m     = 1.435       # m_TO
    E_wr_w     = -0.954      # W/S
    E_wr_A     = 0.2        # aspect ratio A
    E_wr_sweep = -0.702      # wing sweep Lamda
    E_wr_t     = 0.34      # t/c
    E_wr_V     = 0.016        # flight speed V
    E_wr_l     = 0.344       # taper ratio lamda
    E_wr_nz    = 0.686       # nz design positive limit load factor
    E_wr_e     = 0.726       # strut position eta
    
    # Exponents E for strut and juries
    E_wr_m     = 1.556       # m_TO
    E_wr_w     = -1.107      # W/S
    E_wr_A     = 0.885        # aspect ratio A
    E_wr_sweep = -2.516      # wing sweep Lamda
    E_wr_t     = 0.056      # t/c
    E_wr_V     = 0.106        # flight speed V
    E_wr_l     = 1.307       # taper ratio lamda
    E_wr_nz    = 1.148       # nz design positive limit load factor
    E_wr_e     = -4.295       # strut position eta
    E_sj_P     = 46.2        # strut parameter
    
    # Exponents E for aileron efficiency
    E_a_m     = -0.011       # m_TO
    E_a_w     = 0.423      # W/S
    E_a_A     = -0.342        # aspect ratio A
    E_a_sweep = 2.38      # wing sweep Lamda
    E_a_t     = 0.552      # t/c
    E_a_V     = -1.225        # flight speed V
    E_a_l     = -0.075       # taper ratio lamda
    E_a_nz    = 0.522       # nz design positive limit load factor
    E_a_e     = 1.64       # strut position eta
    E_a_P     = -2.634        # strut parameter

    Wred        = 0.15      # weight reduction for the future advanced CFRP
    fold_ratio  = 0.08      # wing weight increment due to wing folding

    sw_chord_ratio  = 0.2       # ratio of strut chord to wing chord at the strut attachment

    n_s             = 0.521798  # strut attachment percentage position

    Pst         = 1-(sw_chord_ratio^0.5*n_s^2)/wing_AR^0.5        # strut parameter

    t_c         = 0.12          # t/c of the wing airfoil

    eta_ail     = C_a*GW^E_a_m*WSR^E_a_w*wing_AR^E_a_A*(m.cos(wing_sweep))^E_a_sweep*(t_c)^E_a_t*




    weight_wing     =                   # wing + strut weight
    #weight_wing     = weight_wing*Units.lb'''



    # RHEA
    weight_wing_SBW     = 5889.4626     #kg
    weight_strut        = 1516.9631     #kg
    #SUGAR real data
    '''weight_wing_SBW     = 7561.3853     #kg
    weight_strut        = 1669.22     #kg'''
    #SUGAR DLR cal
    '''weight_wing_SBW     = 7338.0886     #kg
    weight_strut        = 1697.9703     #kg'''

    weight_wing     = weight_wing_SBW + weight_strut
    #weight_wing     = weight_wing*Units.kg

    return weight_wing
    
   
