## NASA Flops weight methods

#
# Created:  Oct 2019, O. Etlijani
# Modified: Dec 2019, Y. Liu
#           Dec 2019, S. Karpuk
#           Jan 2021, Y. Ma


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
                cab_LE,VARSWP,Machmax,nz,SBWc,Wred_wing,Wred_fold,TFc,mfw,wcr,wct,
                wing_tc_r,wing_tc_t,W_p,fuse_lo,eng_lo,wing_in_sweep,mf_fuel,p_cr,p_SL,
                ETA=[0],chord_length=[0]):
    
    """ Calculate the wing weight of the aircraft based on the FLOPS methods for cantilever wing and 
        semi-empirical methods for Strut braced wing 
    
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
        Machmax        - Max. Mach cruise
        nz             - positive limit load factor  
        SBWc           - strut braced wing configuration
        Wred_fold      - wing folding factor
        TFc            - Twin-fuselage configuration
        mfw            - Estimated wing weight factor (divided by MTOW), is inputs for TF aircraft
        wcr            - wing root chord                                                           m
        wct            - wing tip chord                                                            m
        wing_tc_r      - wing root t/c
        wing_tc_t      - wing tip t/c
        W_p            - weight of propulsion system                                               kg
        fuse_lo        - fuselage spanwise location                                                m
        eng_lo         - engine spanwise location                                                  m
        wing_in_sweep  - inner wing sweep                                                          [radians]
        mf_fuel        - fuel weight factor
        p_cr           - static pressure at cruise altitude
        p_SL           - static pressure sea level

        chord_length   - Local chord length                                                        [m]  # used for detailed wing method, default value 0
             
    Outputs:
      weight_wing                                                                                  [kg]                    
        
    Properties Used:
        N/A
        
    """
    
    # unpack inputs                  
    '''GW           = GW/Units.lb                    # Convert kg to lb
    WSR          = WSR/(Units.lb/Units.ft**2)     # Convert meters squared to ft squared 
    S_wing       = S_wing/Units.ft**2             # Convert from meter squared to ft squared
    wing_span    = wing_span/Units.ft             # Convert meter to ft
    WF           = WF/Units.ft
    OSSPAN       = OSSPAN/Units.ft
    SFLAP        = SFLAP/Units.ft**2
    GLOV         = GLOV/Units.ft**2
    SAFTB        = SAFTB/Units.ft**2
    XLP          = XLP/Units.ft '''

    if SBWc == True:
        # Wing and Strut weight estimation method for Strut-braced wing (SBW) aircraft
        # Source:  P. Chiozzotto G. Initial weight estimate of advanced transport aircraft concepts considering 
        # aeroelastic effects[C]//55th AIAA aerospace sciences meeting. 2017: 0009.

        # Parameters for exponents E and constant C (related to the aircraft configuration and materials)
        # follows parameters valuse are for the SBW configuration with carbon fiber reforced plastic (CFRP)
        
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
        

        # Constant C
        C_c     =   0.00225     # Constant C for wing covers
        C_wr    =   0.209       # Constant C for wing ribs and webs
        C_s    =   0.00101     # Constant C for wing strut and juries
        C_a     =   464         # for aileron efficiency

        # engine relief factor
        ke_c    =   0.99        # covers
        ke_wr   =   0.984       # webs and ribs
        ke_s   =   0.945       # strut and juries
        
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
        E_s_m     = 1.556       # m_TO
        E_s_w     = -1.107      # W/S
        E_s_A     = 0.885        # aspect ratio A
        E_s_sweep = -2.516      # wing sweep Lamda
        E_s_t     = 0.056      # t/c
        E_s_V     = 0.106        # flight speed V
        E_s_l     = 1.307       # taper ratio lamda
        E_s_nz    = 1.148       # nz design positive limit load factor
        E_s_e     = -4.295       # strut position eta
        E_s_P     = 46.2        # strut parameter
        
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

        # ratio of strut chord to wing chord at the strut attachment
        sw_chord_ratio  = 0.1668           # RHEA 0.2; SUGAR 0.1668

        GW = GW * 0.4535924         # TOW should be kg in the follows

        # parameters adjustment
        nz = nz / 1.5
        WSR         = WSR*9.8*(Units.lb/Units.ft**2)     # Convert ft squared to N/m^2
        #Machmax     =   Machmax * 299.278           # equivlent speed!!! covert to m/s
        Machmax = Machmax * 340.29 * m.sqrt(p_cr / p_SL)        # equivlent speed!!! covert to m/s

        AAA = m.cos(wing_sweep)     # not used, only for verify the results

        Pst         = 1-(sw_chord_ratio**0.5*FSTRT**2)/wing_AR**0.5        # strut parameter

        eta_ail     = C_a*GW**E_a_m*WSR**E_a_w*wing_AR**E_a_A*(m.cos(wing_sweep))**E_a_sweep*(wing_TCA)**E_a_t*(
                    Machmax)**E_a_V*(1+wing_TR)**E_a_l*nz**E_a_nz*(1-FSTRT)**E_a_e*(
                        2-FSTRT/(m.cos(wing_sweep))**2)**E_a_P
        
        # calculate wing-box weight penalty factor k_ail
        if eta_ail < 0.5:
            k_ail = (eta_ail/0.5)**-1.1
        else:
            k_ail = 1
        #WSR=4862        ################################################################ only for SUGAR case
        # calculate components group weight
        m_covers = ke_c*C_c*GW**E_c_m*WSR**E_c_w*wing_AR**E_c_A*(m.cos(wing_sweep))**E_c_sweep*(
                    wing_TCA)**E_c_t*Machmax**E_c_V*(1+wing_TR)**E_c_l*nz**E_c_nz*(1-FSTRT)**E_c_e
        
        m_wr = ke_wr*C_wr*GW**E_wr_m*WSR**E_wr_w*wing_AR**E_wr_A*(m.cos(wing_sweep))**E_wr_sweep*(
                wing_TCA)**E_wr_t*Machmax**E_wr_V*(1+wing_TR)**E_wr_l*nz**E_wr_nz*(1-FSTRT)**E_wr_e

        m_strut = ke_s*C_s*GW**E_s_m*WSR**E_s_w*wing_AR**E_s_A*(m.cos(wing_sweep))**E_s_sweep*(
                    wing_TCA)**E_s_t*Machmax**E_s_V*(1+wing_TR)**E_s_l*nz**E_s_nz*(1-FSTRT)**E_s_e*Pst**E_s_P

        # calculate secondary wing structure weight
        m_sec = 0.0443 * GW

        m_wing_before = k_ail*(m_covers+m_wr)+m_sec

        # weight reduction for wing     0.25
        m_wing = m_wing_before * (1-Wred_wing) * (1+Wred_fold)
        # weight reduction for strut    0

        # calculate total weight of wing + strut & juries
        weight_wing     = m_wing + m_strut       # kg
        #weight_wing     = weight_wing*Units.lb

        '''
        # Results comes from Excel file
        # RHEA
        weight_wing_SBW     = 5889.4626     #kg
        weight_strut        = 1516.9631     #kg
        #SUGAR real data
        weight_wing_SBW     = 7561.3853     #kg
        weight_strut        = 1669.22     #kg
        #SUGAR DLR cal
        weight_wing_SBW     = 7338.0886     #kg
        weight_strut        = 1697.9703     #kg

        weight_wing     = weight_wing_SBW + weight_strut
        #weight_wing     = weight_wing*Units.kg'''


    elif TFc == True:
        # Wing weight estimation method for Twin-fuselage aircraft
        # Source: Udin, Wing Mass Formula for Twin Fuselage Aircraft, Journal of Aircraft, 1992
        # A semi-analytical method based on the integrates along the spanwise loads distribution
        # Assumptions: fuel only in the outboard wing; Taper ratio of inner wing is 1. 
        # All units are IS standard: m, kg, N, ...
        
        #inputs
        sf_ail          = 0.04          # ailerons area factor (divided by wing area)
        sf_flaps        = 0.18          # flaps area factor (divided by wing area)
        E_T             = 1.3           # effective airfoil thickness coefficient
        g               = 9.8           # gravitational accelerate, N/m^2

        # flaps factors
        k_lef           = 3.5           # Krueger flap is 2.5, others is 3.5
        k_sup           = 1.6           # flap motion support, simple hinge external support 1.0; link/track end support 1.2; Fowler flaps with hooked track 1.6
        k_slot          = 1.5           # the number of slots, 1.0 single-slotted, 1.5 double-slot, 2.0 double-slotted flaps with articulating vanes
        Ome_ref         = 56            # N/m^2, reference specific weight, is defined equal to the hypothetical weight of a thin section built up from two Al-alloy skins with one millimetre thickness each.
        W_ref           = 1000000       # N

        # Materials properties
        rho             = 1550.0746     # kg/m^3, metarial density of wing structures, CFPR is 1550, comes from D8 report; 7075 Al alloy is 2800
        Sig_u           = 410927517.2   # Pa, N/m^2, ultimate direct stress, CFPR is 410927517.2, comes from D8 report; 7075 Al is 310251000
        Sig_us          = 441264448     # Pa, N/m^2, ultimate shear stress, CFPR is 441264448, comes from D8 report; 7075 Al is 310251000

        # internal parameter calculations
        nz = nz / 1.5
        W_p = W_p * 1.12                # 1.12 is for the nacelle weight
        half_span  = wing_span / 2
        T_r = wcr * wing_tc_r           # absolute wing root thickness
        T_t = wct * wing_tc_t           # absolute wing tip thickness
        wing_sweep = wing_sweep * 0.92  # outboard wing sweep, transfer 1/4c wing sweep to 1/2c wing sweep

        mf_eng = W_p / GW               # engine weight factor
        z_f = fuse_lo / half_span       # relative fuselage location
        z_e = eng_lo / half_span        # relative engine location

        h   = T_t / T_r                 # outboard wing thickness taper ratio
        PSI = wing_TR ** 2 *(wing_tc_t / wing_tc_r)         # outboard wing airfoil area taper ratio

        # Factors for wing-box weight
        k_sl            = 1.2           # service life factor, ultimate stress divided by panel fatigue stress
        k_man           = 1.4           # manufacturing factor, for conventional aircraft, 1.62-2.05; for advanced technology aircraft, 1.3
        k_tw            = 1 + m.sqrt(wing_AR*(2*z_f+(1-z_f)*(1+wing_TR))) * (0.0225/m.cos(wing_in_sweep)*m.sqrt(z_f/2) + 
                           0.015*(1+2*wing_TR)*m.sqrt(1-z_f)/((wing_TR+1)**1.5 * m.cos(wing_sweep)))
        
        # Mass required by shear force
        K_Qi = 1/m.cos(wing_in_sweep) * (z_f**2 / (z_f*(1-wing_TR)+wing_TR+1) - 
                z_f**2*(2+1.2*wing_TR)*mfw/(z_f*(1-wing_TR)+4*wing_TR+4))

        K_Qo = 1/m.cos(wing_sweep) * ((1-z_f)**2 /3 * ((1+2*wing_TR)/(z_f*(1-wing_TR)+wing_TR+1) -
                (1+2*PSI)*mf_fuel/(PSI+1)/(1-z_f)-(4-5.6*wing_TR)*mfw/(z_f*(1-wing_TR)+4*wing_TR+4)) -
                (z_e-z_f)*mf_eng)
        
        m_Q  = rho*g*nz/2/Sig_us * m.sqrt(GW*wing_AR/WSR) * (K_Qi + K_Qo)

        # Mass required by bending moment
        K_Mi = z_f/m.cos(wing_in_sweep) * (((1-z_f)**2*(1+2*wing_TR)-2*z_f**2)/3/(z_f*(1-wing_TR)+wing_TR+1) -
                (1-z_f)*(1+2*PSI)*mf_fuel/(3*PSI+3) - 
                ((4/3+1.87*wing_TR)*(1-z_f)**2-17/12*z_f**2*(1+0.6*wing_TR))*mfw / (z_f*(1-wing_TR)+4*(wing_TR+1)) -
                z_e * mf_eng)
        
        A1 = (1-z_f)**2/3/(1-h)**3 * (1/3 - 3/2*h + 3*h**2 - 11/6*h**3 + h**3*m.log10(h))
        A2 = (1-wing_TR)/(z_f*(1-wing_TR)+wing_TR+1) - (1-PSI)*mf_fuel/(1+PSI)/(1-z_f) - (4+0.8*wing_TR)*mfw/(z_f*(1-wing_TR)+4+4*wing_TR)
        A3 = ((1-z_f)/(1-h))**2 * (1/2-2*h+3/2*h**2-h**2*m.log10(h))
        A4 = wing_TR/(z_f*(1-wing_TR)+wing_TR+1) - PSI*mf_fuel/(PSI+1)/(1-z_f) - 1.6*wing_TR*mfw/(z_f*(1-wing_TR)+4*wing_TR+4)
        A5 = mf_eng * (z_e-z_f-(h*(1-z_f)/(h-1)-1+z_e)*m.log10((1-z_e)/(1-z_f)*(1-h)+h))
        K_Mo = (1-z_f)/(1-h)/m.cos(wing_sweep) * (A1 * A2 + A3 * A4 - A5)

        m_M = rho*g*nz/2/Sig_u * wing_AR**1.5 * m.sqrt(GW/WSR) * E_T * (wing_TR+1)*(1-z_f)/2/wing_tc_r * (K_Mi + K_Mo)

        # wing-box weight factor: bending moment + shear force
        m_wingbox   =   k_sl * k_man * k_tw *(m_Q + m_M)

        # wing secondary components weight factor
        m_rib   =   0.15 * m_M
        m_ail   =   0.03 * sf_ail
        m_sk    =   3 / WSR
        
        m_lef   =   k_lef / WSR
            # Only trailing edge flaps weight method from Torenbeek 2013 book
            # *0.8 means the weight reduction from the Al to CFPR, the original method is for Al materials
        M_tef   =   1.7 * k_sup * k_slot * Ome_ref * (1 + (GW * 9.8 / W_ref)) * S_wing * sf_flaps / g * 0.8
        m_tef   =   M_tef / GW

        m_of    =   30 * sf_flaps / WSR
        m_flap  =   m_lef + m_tef + m_of

        m_wingsec   =   m_rib + m_ail + m_sk + m_flap

        # Total wing factor
        mf_w    =   m_wingbox + m_wingsec
        mf_w    =   mf_w * (1-Wred_wing) * (1+Wred_fold)
        
        # Total wing weight 
        weight_wing     =   mf_w * GW               # kg

        




    else:
        # Cantilever wing weight
        #Calculate weight of wing using FLOPS methods

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
        weight_wing     = weight_wing * (1-Wred_wing) * (1+Wred_fold)
        weight_wing     = weight_wing*Units.lb
        




    return weight_wing
    
   
