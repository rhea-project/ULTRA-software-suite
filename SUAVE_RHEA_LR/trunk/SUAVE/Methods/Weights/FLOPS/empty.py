## @ingroup Methods-Weights-FLOPS
# empty.py
# 
# Created:  Dec 2019, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import numpy as np
import SUAVE
from SUAVE.Core                                                                         import Units, Data
from SUAVE.Methods.Weights.FLOPS.systems                                                import systems                   
from SUAVE.Methods.Weights.FLOPS.weight_payload_items                                   import weight_payload                   as payload
from SUAVE.Methods.Weights.FLOPS.weight_wing                                            import weight_wing                      as weight_wing
from SUAVE.Methods.Weights.FLOPS.weight_horizontal_tail                                 import weight_ht                        as weight_horizontal_tail
from SUAVE.Methods.Weights.FLOPS.weight_vertical_tail                                   import weight_vt                        as weight_vertical_tail
from SUAVE.Methods.Weights.FLOPS.weight_fuselage                                        import weight_fuselage                  as weight_fuselage
from SUAVE.Methods.Weights.FLOPS.weight_landing_gear                                    import weight_lg                        as weight_lg
from SUAVE.Methods.Weights.FLOPS.weight_fin                                             import weight_fin                       as weight_fin
from SUAVE.Methods.Weights.FLOPS.weight_canard                                          import weight_canard                    as weight_canard
from SUAVE.Methods.Weights.FLOPS.weight_paint                                           import weight_paint                     as weight_paint  
from SUAVE.Methods.Weights.FLOPS.Propulsion_systems.engine_weight                       import weight_engine                    as weight_engine 
from SUAVE.Methods.Weights.FLOPS.Propulsion_systems.miscellaneous_propulsion_systems    import weight_engine_miscellaneous      as weight_engine_miscellaneous
from SUAVE.Methods.Weights.FLOPS.Propulsion_systems.thrust_reversers                    import weight_thrust_reversers          as weight_thrust_reversers
from SUAVE.Methods.Weights.FLOPS.weight_nacelles_air_induciton                          import weight_nac                       as weight_nac
from SUAVE.Methods.Weights.Correlations.Common                                          import wing_main                        as wing_main
from SUAVE.Methods.Weights.EMWET                                                        import wing_main                        as EMWET_wing

import warnings

# ----------------------------------------------------------------------
#  Empty
# ----------------------------------------------------------------------

## @ingroup Methods-Weights-FLOPS
def empty(vehicle):
    """ This scripts calculates the empty weight based on the FLOPS Method.        

    Assumptions:
        calculated aircraft weight from the FLOPS method
    
    Source:
        N/A
        
    Inputs:
      engine - a data dictionary with the fields:                    
          thrust_sls - sea level static thrust of a single engine                [Newtons]

      wing - a data dictionary with the fields:
          gross_area - wing gross area                                           [meters**2]
          span - span of the wing                                                [meters]
          taper - taper ratio of the wing                                        [dimensionless]
          t_c - thickness-to-chord ratio of the wing                             [dimensionless]
          sweep - sweep angle of the wing                                        [radians]
          mac - mean aerodynamic chord of the wing                               [meters]
          r_c - wing root chord                                                  [meters]

      aircraft - a data dictionary with the fields:                    
          Nult - ultimate load of the aircraft                                   [dimensionless]
          Nlim - limit load factor at zero fuel weight of the aircraft           [dimensionless]
          TOW - maximum takeoff weight of the aircraft                           [kilograms]
          zfw - maximum zero fuel weight of the aircraft                         [kilograms]
          num_eng - number of engines on the aircraft                            [dimensionless]
          num_pax - number of passengers on the aircraft                         [dimensionless]
          wt_cargo - weight of the bulk cargo being carried on the aircraft      [kilograms]
          num_seats - number of seats installed on the aircraft                  [dimensionless]
          ctrl - specifies if the control system is "fully powered", "partially powered", or not powered [dimensionless]
          ac - determines type of instruments, electronics, and operating items based on types: 
              "short-range", "medium-range", "long-range", "business", "cargo", "commuter", "sst"        [dimensionless]
          w2h - tail length (distance from the airplane c.g. to the horizontal tail aerodynamic center)  [meters]

      fuselage - a data dictionary with the fields:
          area - fuselage wetted area                                            [meters**2]
          diff_p - Maximum fuselage pressure differential                        [Pascal]
          width - width of the fuselage                                          [meters]
          height - height of the fuselage                                        [meters]
          length - length of the fuselage                                        [meters]                     

      horizontal
          area - area of the horizontal tail                                     [meters**2]
          span - span of the horizontal tail                                     [meters]
          sweep - sweep of the horizontal tail                                   [radians]
          mac - mean aerodynamic chord of the horizontal tail                    [meters]
          t_c - thickness-to-chord ratio of the horizontal tail                  [dimensionless]
          exposed - exposed area ratio for the horizontal tail                   [dimensionless]

      vertical
          area - area of the vertical tail                                       [meters**2]
          span - sweight = weight * Units.lbpan of the vertical                  [meters]
          t_c - thickness-to-chord ratio of the vertical tail                    [dimensionless]
          sweep - sweep angle of the vertical tail                               [radians]
          t_tail - factor to determine if aircraft has a t-tail, "yes"           [dimensionless]
          
      settings.weight_reduction_factors.
          main_wing                                                              [dimensionless] (.1 is a 10% weight reduction)
          empennage                                                              [dimensionless] (.1 is a 10% weight reduction)
          fuselage                                                               [dimensionless] (.1 is a 10% weight reduction)
          

    Outputs:
        output - a data dictionary with fields:
            wt_payload - weight of the passengers plus baggage and paid cargo    [kilograms]
            wt_pax - weight of all the passengers                                [kilogram]
            wt_bag - weight of all the baggage                                   [kilogram]
            wt_fuel - weight of the fuel carried                                 [kilogram]
            wt_empty - operating empty weight of the aircraft                    [kilograms]
 
    Properties Used:
        N/A
    """     

    # Unpack inputs
    settings   = vehicle.settings
    p_h        = vehicle.cruise_pressure
    p_0        = vehicle.SL_pressure
    Mc         = vehicle.cruise_mach
    Mmax       = vehicle.max_cruise_mach
    cr_range   = vehicle.range 
    n_z        = vehicle.envelope.ultimate_load
    if Mmax < 1:
        cruise_type = 1
    else:
        cruise_type = 2
            
    Nult       = vehicle.envelope.ultimate_load
    Nlim       = vehicle.envelope.limit_load
    TOW        = vehicle.mass_properties.max_takeoff
    LW         = vehicle.mass_properties.max_landing
    wt_zf      = vehicle.mass_properties.max_zero_fuel
    num_crew   = vehicle.crew
    num_pax    = vehicle.passengers
    num_first  = vehicle.passengers_first
    num_biz    = vehicle.passengers_business
    num_econ   = vehicle.passengers_economy
    wt_cargo   = vehicle.mass_properties.cargo
    wt_cargo_w = vehicle.mass_properties.wing_cargo           
    wt_cargo_f = vehicle.mass_properties.fuselage_cargo
    FULAUX     = vehicle.mass_properties.auxilary_fuel
    ac_type    = vehicle.systems.accessories    
    
    # Propulsion Group ----------------------------------------------------------------------------------------    
    propulsor_name = list(vehicle.propulsors.keys())[0]                     #obtain the key for the propulsor for assignment purposes
    
    propulsors      = vehicle.propulsors[propulsor_name]
    num_eng         = propulsors.number_of_engines
    engine_pos      = propulsors.engine_position
    engine_Yloc     = propulsors.origin[0][1]
    engine_d        = propulsors.nacelle_diameter
    Swet_eng        = propulsors.areas.wetted 
    eng_w           = propulsors.wing_mounted_engines
    eng_f           = propulsors.fuselage_mounted_engines
    eng_origin      = propulsors.origin 
    eng_weight_unit = np.zeros_like(eng_origin[0][:])                     # initialize the engine weight matrix for each nagine 
                                                                          # (starts with wing-mounted engines and ends witn fuselage/moundet)

    # Propulsion and all related Groups ---------------------------------------------------------------------------------------    

    if propulsor_name=='turbofan' or propulsor_name=='Turbofan':
        # thrust_sls should be sea level static thrust.       
        thrust_sls   = propulsors.thrust.total_design
        thrust_base  = propulsors.thrust.baseline
        WENGB        = propulsors.base_engine_weight
        WINLB        = propulsors.inlet_weight
        WNOZB        = propulsors.nozzle_weight
        WPMISC       = propulsors.miscellaneous_weight
        EEXP         = propulsors.weight_scaling
        inl_inc      = propulsors.nozzle_included 
        noz_inc      = propulsors.inlet_included
        l_eng        = propulsors.engine_length
        d_eng        = propulsors.nacelle_diameter

        # Engine
        # Note: the baseline engine thrust is assumed as static sea-level thrust from the propulsion
        #       definition. If otherwise, then specify it in the vehicle set-up file
        if thrust_base == 0:
            thrust_base = thrust_sls

        wt_engine = weight_engine(vehicle.type,eng_w,eng_f,thrust_sls/num_eng,thrust_base/num_eng, \
                                  WENGB,EEXP,inl_inc,noz_inc,1,1,WINLB,WNOZB)
        # for future engine weight reduction, calculated refer to SUGAR
        #wt_engine = wt_engine * 0.85

        # assign engine weights to the wing or the fuselage
        wt_engine_w = wt_engine / num_eng 
        wt_engine_f = wt_engine / num_eng

        cond_ind = 0
        if eng_w != 0:
            eng_weight_unit[0:eng_w] = wt_engine_w
            cond_ind = eng_w
        if eng_f != 0:
            eng_weight_unit[cond_ind:cond_ind+eng_f] = wt_engine_f 
            cond_ind += eng_f

        # Miscelaneous
        wt_eng_misc = weight_engine_miscellaneous(eng_w,eng_f,vehicle.type,thrust_sls/num_eng,engine_d, \
                                                                WPMISC,num_crew,Mmax)

        # Thrust_reversers
        if propulsors.thrust_reverser== True:
            wt_thrust_rev = weight_thrust_reversers(eng_w,eng_f,thrust_sls/num_eng)
        else:
            wt_thrust_rev = 0

    elif propulsor_name=='turbofan_hybrid' or propulsor_name=='Turbofan_hybrid': 
        # thrust_sls should be sea level static thrust.
        thrust_sls        = propulsors.thrust.total_design 
        thrust_base       = propulsors.thrust.baseline
        thrust_sls_jet    = propulsors.thrust.total_design * (1-propulsors.thrust.hybridization)
        thrust_sls_el     = propulsors.thrust.total_design * propulsors.thrust.hybridization
        thrust_base_jet   = propulsors.thrust.baseline     * (1-propulsors.thrust.hybridization)
        thrust_base_el    = propulsors.thrust.baseline * propulsors.thrust.hybridization
        WENGB             = propulsors.base_engine_weight
        WINLB             = propulsors.inlet_weight
        WNOZB             = propulsors.nozzle_weight
        WPMISC            = propulsors.miscellaneous_weight
        EEXP              = propulsors.weight_scaling
        noz_inc           = propulsors.nozzle_included 
        inl_inc           = propulsors.inlet_included
        num_jet_eng       = propulsors.jet_engines
        num_jet_eng_w     = propulsors.wing_mounted_jet_engines
        num_jet_eng_f     = propulsors.fuselage_mounted_jet_engines 
        num_el_eng        = propulsors.electric_engines
        num_el_eng_w      = propulsors.wing_mounted_electric_engines 
        num_el_eng_f      = propulsors.fuselage_mounted_electric_engines
        l_eng_jet         = propulsors.engine_length
        d_eng_jet         = propulsors.nacelle_diameter
        l_eng_el          = propulsors.engine_length
        d_eng_el          = propulsors.nacelle_diameter

        # assign positions in the propulsion weights matrix for each type of engines


        # Correction for extreme cases when the aircraft is only jet-powered or electrically-powered
        if propulsors.thrust.hybridization == 0:
            num_jet_eng   = num_eng;                         num_el_eng   = 0                            
            num_jet_eng_w = num_jet_eng_w + num_el_eng_w;    num_el_eng_w = 0      
            num_jet_eng_f = num_jet_eng_f + num_el_eng_f;    num_el_eng_f = 0   
        elif propulsors.thrust.hybridization == 1:
            num_el_eng   = num_eng;                         num_jet_eng   = 0                            
            num_el_eng_w = num_jet_eng_w + num_el_eng_w;    num_jet_eng_w = 0      
            num_el_eng_f = num_jet_eng_f + num_el_eng_f;    num_jet_eng_f = 0               

        # Gas turbine Engine
        if thrust_base_jet == 0:
            thrust_base = thrust_sls

        wt_engine1 = weight_engine(vehicle.type,num_jet_eng_w,num_jet_eng_f,thrust_sls_jet/num_jet_eng,thrust_base_jet/num_jet_eng, \
                                   WENGB,EEXP,inl_inc,noz_inc,1,1,WINLB,WNOZB)

        wt_engine_w1 = wt_engine1 / num_jet_eng
        wt_engine_f1 = wt_engine1 / num_jet_eng 

        # Miscelaneous
        wt_eng_misc1 = weight_engine_miscellaneous(num_jet_eng_w,num_jet_eng_f,vehicle.type,thrust_sls_jet/num_jet_eng,engine_d, \
                                                                WPMISC,num_crew,Mmax)

        # Thrust_reversers
        if propulsors.thrust_reverser== True:
            wt_thrust_rev1 = weight_thrust_reversers(num_jet_eng_w,num_jet_eng_f,thrust_sls/num_jet_eng)
        else:
            wt_thrust_rev1 = 0

        # Electric fan engine
        T_W          = propulsors.motor.power_to_weight       
        power_sls    = propulsors.max_power_electric/1000   # Sea-level electric engine power in kW
        wt_eng_misc  = propulsors.miscellaneous_weight 

        wt_engine2 = 2*power_sls / T_W          # Assume installation weight factor of the motor weight of 2

        wt_engine_w2 = wt_engine2 / num_el_eng
        wt_engine_f2 = wt_engine2 / num_el_eng

        # Thrust_reversers
        wt_thrust_rev2 = 0

        # Miscelaneous
        wt_eng_misc2 = wt_eng_misc

        # Add gas turbine and electric engines
        wt_engine     = wt_engine1     + wt_engine2
        wt_eng_misc   = wt_eng_misc1   + wt_eng_misc2
        wt_thrust_rev = wt_thrust_rev1 + wt_thrust_rev2

        # Assign individual engine weights to each position
        cond_ind = 0
        if num_jet_eng_w != 0:
            eng_weight_unit[0:num_jet_eng_w] = wt_engine_w1
            cond_ind = num_jet_eng_w
        if num_el_eng_w != 0:
            eng_weight_unit[cond_ind:cond_ind+num_el_eng_w] = wt_engine_w2
            cond_ind += num_el_eng_w
        if num_jet_eng_f != 0:
            eng_weight_unit[cond_ind:cond_ind+num_jet_eng_f] = wt_engine_f1 
            cond_ind += num_jet_eng_f
        if num_el_eng_f != 0:
            eng_weight_unit[cond_ind:cond_ind+num_el_eng_f] = wt_engine_f2 
            cond_ind += num_el_eng_f

    else:
        # Alternate engine data definition as in FLOPS. Mostly applicable to electric aircraft
        # Engine
        T_W          = propulsors.motor.thrust_to_weight        
        thrust_sls   = propulsors.thrust.total_design
        wt_eng_misc  = propulsors.miscellaneous_weight 

        wt_engine = 2*thrust_sls / T_W          # Assume installation weight of the motor weight

        # Thrust_reversers
        wt_thrust_rev = 0

        # Miscelaneous
        wt_eng_misc = wt_eng_misc

    wt_propulsion = wt_engine + wt_eng_misc + wt_thrust_rev
    propulsors.mass_properties.mass  = wt_propulsion
        

    if wt_propulsion==0:
        warnings.warn("Propulsion mass= 0 ;e there is no Engine Weight being added to the Configuration", stacklevel=1)
        
    # -------------------------------------------------------------------------------------------------------    

    # Fuselage group ----------------------------------------------------------------------------------------      
    if 'fuselage' not in vehicle.fuselages:        
        wt_fuselage = 0
        S_fus       = 0
    else:
        Wfus        = vehicle.fuselages['fuselage'].reduction_factor
        S_fus       = vehicle.fuselages['fuselage'].areas.wetted
        cabin_area  = vehicle.fuselages['fuselage'].areas.cabin_floor
        diff_p_fus  = vehicle.fuselages['fuselage'].differential_pressure
        w_fus       = vehicle.fuselages['fuselage'].width
        h_fus       = vehicle.fuselages['fuselage'].heights.maximum
        l_fus       = vehicle.fuselages['fuselage'].lengths.total
        l_cab       = vehicle.fuselages['fuselage'].lengths.cabin
        n_fus       = vehicle.fuselages['fuselage'].number_of_fuselages
        num_seats   = vehicle.fuselages['fuselage'].number_coach_seats
        CARGF       = vehicle.fuselages['fuselage'].CARGF
        CAYF        = vehicle.fuselages['fuselage'].fuselage_number_factor
        var_sweep   = vehicle.wings['main_wing'].sweeps.variable
        RSPCHD      = vehicle.wings['main_wing'].aft_body.aft_spar_location_center
        sweep_w     = vehicle.wings['main_wing'].sweeps.quarter_chord

        if n_fus != 1:
            CAYF = 0.5
        else:
            CAYF = 1
            
        wt_fuselage = weight_fuselage(l_fus,l_cab,sweep_w,eng_f,CARGF,n_fus,w_fus,h_fus,TOW, \
                                      vehicle.type,cabin_area,var_sweep,Nult,RSPCHD,p_h,p_0,Mc) 
        wt_fuselage = wt_fuselage*(1.-Wfus)

        vehicle.fuselages['fuselage'].mass_properties.mass = wt_fuselage
    # -------------------------------------------------------------------------------------------------------    

    # Wing Group --------------------------------------------------------------------------------------------            
    S_gross_w  = vehicle.reference_area
    if 'main_wing' not in vehicle.wings:
        wt_wing  = 0.0
        wing_c_r = 0.0
        Swet_w   = 0.0
        warnings.warn("There is no Wing Weight being added to the Configuration", stacklevel=1)
        
    else:
        Wred       = vehicle.wings['main_wing'].reduction_factor
        fold_ratio = vehicle.wings['main_wing'].folding_penalty
        AR         = vehicle.wings['main_wing'].aspect_ratio
        Swet_w     = vehicle.wings['main_wing'].areas.wetted
        SFLAP      = vehicle.wings['main_wing'].areas.flap  
        WSR        = vehicle.wings['main_wing'].wing_loading
        GLOV       = vehicle.wings['main_wing'].GLOV 
        b          = vehicle.wings['main_wing'].spans.projected
        lambda_w   = vehicle.wings['main_wing'].taper
        delta_w    = vehicle.wings['main_wing'].dihedral
        t_c_w      = vehicle.wings['main_wing'].thickness_to_chord
        sweep_w    = vehicle.wings['main_wing'].sweeps.quarter_chord
        var_sweep  = vehicle.wings['main_wing'].sweeps.variable
        mac_w      = vehicle.wings['main_wing'].chords.mean_aerodynamic
        wing_c_r   = vehicle.wings['main_wing'].chords.root
        wing_c_t   = vehicle.wings['main_wing'].chords.tip
        comp_f     = vehicle.wings['main_wing'].composite_fraction
        S_aftb     = vehicle.wings['main_wing'].aft_centerbody_area 
        OSSPAN     = vehicle.wings['main_wing'].outerwing.spans
        FSTRT      = vehicle.wings['main_wing'].FSTRT
        FAERT      = vehicle.wings['main_wing'].FAERT
        PCTL       = vehicle.wings['main_wing'].PCTL
        RSPCHD     = vehicle.wings['main_wing'].aft_body.aft_spar_location_center
        RSPSOB     = vehicle.wings['main_wing'].aft_body.aft_spar_location_side
        OSSPAN     = vehicle.wings['main_wing'].outerwing.spans
        XLP        = vehicle.wings['main_wing'].cabin_length
        cab_LE     = vehicle.wings['main_wing'].cabin_sweep_leading_edge
        SBWcon     = vehicle.wings['main_wing'].SBW
        TFcon      = vehicle.wings['main_wing'].TF
        mf_w       = vehicle.wings['main_wing'].mf_wing
        t_c_r      = vehicle.wings['main_wing'].t_c_ratio_root          # root t/c
        t_c_t      = vehicle.wings['main_wing'].t_c_ratio_tip          # tip t/c
        fuse_location = vehicle.wings['main_wing'].fuselage_location
        eng_location  = vehicle.wings['main_wing'].engine_location
        sweep_in_w    = vehicle.wings['main_wing'].sweeps_inner_half_chord
        wf_fuel    = vehicle.wings['main_wing'].fuel_factor


    # Unpack and run EMWET if the flag is True
        if vehicle.wings['main_wing'].EMWET is True:
            wing_eng_num = propulsors.wing_mounted_engines                              # Assume each side has an equal amount of engines        
            wt_wing_mat = EMWET_wing(vehicle,eng_origin,wing_eng_num,eng_weight_unit)
            wt_wing = wt_wing_mat[0]
        else:
            wt_wing = weight_wing(1,TOW,WSR,S_gross_w,sweep_w,lambda_w,t_c_w,b,AR,GLOV, \
                                settings.wing_weight_method,settings.wing_weight_method, \
                                vehicle.type,w_fus,OSSPAN,FSTRT,FAERT,PCTL,SFLAP,eng_w,comp_f, \
                                eng_f,S_aftb,RSPCHD,RSPSOB,XLP,Nult,CAYF,cab_LE,var_sweep,Mmax,n_z, \
                                SBWcon,Wred,fold_ratio,TFcon,mf_w,wing_c_r,wing_c_t,t_c_r,t_c_t, \
                                wt_propulsion,fuse_location,eng_location,sweep_in_w,wf_fuel,p_h,p_0)
        #wt_wing    = wing_main.wing_main(S_gross_w,b,lambda_w,t_c_w,sweep_w,Nult,TOW,wt_zf)

        #wt_wing    = wt_wing*(1.-Wred)*(1+fold_ratio)   
        #wt_wing    = wt_wing*(1+fold_ratio)
        vehicle.wings['main_wing'].mass_properties.mass = wt_wing            
    # -------------------------------------------------------------------------------------------------------    

    # Tail group------------------------------------------------------------------------------------------------------    
    # Horizontal tail
    if 'horizontal_stabilizer' not in vehicle.wings:
        wt_tail_horizontal = 0.0
        S_h = 0.0
        warnings.warn("There is no Horizontal Tail Weight being added to the Configuration", stacklevel=1)
        
    else:    
        Wred_h    = vehicle.wings['horizontal_stabilizer'].reduction_factor
        S_h       = vehicle.wings['horizontal_stabilizer'].areas.reference
        Swet_h    = vehicle.wings['horizontal_stabilizer'].areas.wetted
        lambda_h  = vehicle.wings['horizontal_stabilizer'].taper
        
        wt_tail_horizontal = weight_horizontal_tail(S_h,TOW,lambda_h,Nult,Mc,p_h,p_0,vehicle.type)   
        wt_tail_horizontal = wt_tail_horizontal*(1.-Wred_h)
        vehicle.wings['horizontal_stabilizer'].mass_properties.mass = wt_tail_horizontal

    # Vertical tail
    if 'vertical_stabilizer' not in vehicle.wings:   
        output_3                  = Data()
        output_3.wt_tail_vertical = 0.0
        output_3.wt_rudder        = 0.0
        S_v                       = 0.0
        Swet_v                    = 0.0
        warnings.warn("There is no Vertical Tail Weight being added to the Configuration", stacklevel=1)            
    else:    
        Wred_v       = vehicle.wings['vertical_stabilizer'].reduction_factor 
        S_v          = vehicle.wings['vertical_stabilizer'].areas.reference
        Swet_v       = vehicle.wings['vertical_stabilizer'].areas.wetted
        lambda_v     = vehicle.wings['vertical_stabilizer'].taper
        AR_v         = vehicle.wings['vertical_stabilizer'].aspect_ratio
        b_v          = vehicle.wings['vertical_stabilizer'].spans.projected
        t_c_v        = vehicle.wings['vertical_stabilizer'].thickness_to_chord
        sweep_v      = vehicle.wings['vertical_stabilizer'].sweeps.quarter_chord
        NVERT        = vehicle.wings['vertical_stabilizer'].number_of_wings  
        HHT          = vehicle.wings['vertical_stabilizer'].wing_mount_coef           
        
        wt_tail_vertical = weight_vertical_tail(TOW,lambda_v,NVERT,S_v,Nult,AR_v,sweep_v,HHT,t_c_v,p_h,p_0,Mc,vehicle.type)
        wt_tail_vertical = wt_tail_vertical*(1.-Wred_v)
        
        vehicle.wings['vertical_stabilizer'].mass_properties.mass = wt_tail_vertical
    # -------------------------------------------------------------------------------------------------------    

    # Fin ---------------------------------------------------------------------------------------------------
    if 'fin' not in vehicle:
        Swet_fin   = 0
        wt_fin_tot = 0
    else:
        S_fin      = vehicle.wings['fin'].areas.reference
        Swet_fin   = vehicle.wings['fin'].areas.wetted
        lambda_f   = vehicle.wings['fin'].taper
        NFIN       = vehicle.wings['fin'].number_of_tails
        wt_fin_tot = weight_fin(TOW,S_fin,lambda_f,NFIN)
        wt_fin_tot = wt_fin_tot*(1.-wt_factors.empennage)        
    # -------------------------------------------------------------------------------------------------------

    # Canard ------------------------------------------------------------------------------------------------
    if 'canard' not in vehicle:
        Swet_can       = 0
        wt_cannard_tot = 0
    else:
        S_can         = vehicle.wings['canard'].areas.reference
        Swet_can      = vehicle.wings['canard'].areas.wetted
        lambda_c      = vehicle.wings['canard'].taper
        wt_canard_tot = weight_canard(TOW,S_can,lambda_c)
        wt_canard_tot = wt_canard_tot*(1.-wt_factors.empennage)         
    # -------------------------------------------------------------------------------------------------------
    
    # Landing gear group -------------------------------------------------------------------------------------           
    length_str      = vehicle.landing_gear.main_strut_length
    length_nose_str = vehicle.landing_gear.nose_strut_length
    
    wt_landing_gear, wt_nose_gear, wt_main_gear = weight_lg(vehicle.type,engine_pos,delta_w,engine_Yloc,LW,engine_d,w_fus,l_fus, \
                                                            cr_range,length_nose_str,length_str,cruise_type,TOW)  
    # -------------------------------------------------------------------------------------------------------       

    # Systems group -----------------------------------------------------------------------------------------   
    output_2 = systems(vehicle) 
    # -------------------------------------------------------------------------------------------------------    
   
    # Paint -------------------------------------------------------------------------------------------------    
    paint_dens = vehicle.paint.area_density
    wt_paint   = weight_paint(paint_dens,Swet_w,Swet_h,Swet_v,S_fus,Swet_eng,Swet_can,Swet_fin,n_fus,NVERT)
    # -------------------------------------------------------------------------------------------------------    

    # Battery group -----------------------------------------------------------------------------------------
    if propulsors.type == 'electric' or propulsors.type == 'hybrid':
        wt_battery      = propulsors.battery.mass_properties.mass
        wt_bat_payload  = 0 #propulsors.payload.mass_properties.mass
    else:
        wt_battery     = 0
        wt_bat_payload = 0
    # -------------------------------------------------------------------------------------------------------

    # Electric cables -----------------------------------------------------------------------------------------
    if propulsors.type == 'electric':
        wt_cables = propulsors.cables.density * propulsors.cables.total_length 
    elif propulsors.type == 'hybrid':
        if propulsors.thrust.hybridization != 0:
            wt_cables = propulsors.cables.density * propulsors.cables.total_length 
        else:
            wt_cables = 0
    else:
        wt_cables = 0

    # Nacelle air inuction Group ---------------------------------------------------------------------------
    if propulsors.type == 'hybrid':
        if propulsors.thrust.hybridization == 0:
            wt_air_ind = weight_nac(vehicle.type,num_eng,thrust_sls/num_eng,thrust_base/num_eng,eng_f,l_fus,h_fus,l_eng_jet,d_eng_jet)
        elif propulsors.thrust.hybridization == 1:
            wt_air_ind = weight_nac(vehicle.type,num_eng,thrust_sls/num_eng,thrust_base/num_eng,eng_f,l_fus,h_fus,l_eng_el,d_eng_el) 
        else:
            wt_air_ind1 = weight_nac(vehicle.type,num_jet_eng,thrust_sls_jet/num_jet_eng,thrust_base_jet/num_jet_eng,num_jet_eng_f,l_fus,h_fus,l_eng_jet,d_eng_jet)
            wt_air_ind2 = weight_nac(vehicle.type,num_el_eng,thrust_sls_el/num_el_eng,thrust_base_el/num_el_eng,num_el_eng_f,l_fus,h_fus,l_eng_el,d_eng_el)
            wt_air_ind  = wt_air_ind1+wt_air_ind2
    elif propulsors.type != 'electric':
        wt_air_ind = weight_nac(vehicle.type,num_eng,thrust_sls/num_eng,thrust_base/num_eng,eng_f,l_fus,h_fus,l_eng,d_eng)
    else:
        wt_air_ind = 0
    #--------------------------------------------------------------------------------------------------------
    
    # Calculate the equipment empty weight of the aircraft
    wt_empty = (wt_wing + wt_fuselage + wt_tail_horizontal + wt_tail_vertical + wt_fin_tot + wt_cannard_tot + \
                wt_landing_gear + output_2.wt_systems + wt_paint + wt_propulsion + wt_battery + wt_air_ind + wt_cables)

    
    # packup outputs
    output = Data()
    output.empty             = wt_empty
    output.payload           = payload(vehicle.type,cr_range,num_first,num_biz,num_econ,wt_cargo_w,wt_cargo_f,FULAUX)
    output.wing              = wt_wing
    output.fuselage          = wt_fuselage
    output.propulsion        = wt_propulsion
    output.nacelles          = wt_air_ind
    output.landing_gear      = wt_landing_gear
    output.horizontal_tail   = wt_tail_horizontal
    output.vertical_tail     = wt_tail_vertical
    output.battery           = wt_battery
    output.battery_payload   = wt_bat_payload
    output.cables            = wt_cables
    output.paint             = wt_paint
    output.systems           = output_2.wt_systems       
    output.systems_breakdown = Data()
    output.systems_breakdown.apu               = output_2.wt_apu   
    output.systems_breakdown.avionics          = output_2.wt_avionics        
    output.systems_breakdown.electrical        = output_2.wt_electrical    
    output.systems_breakdown.furnish           = output_2.wt_furnish
    output.systems_breakdown.hydraulics        = output_2.wt_hyd_pnu
    output.systems_breakdown.instruments       = output_2.wt_instruments
    output.systems_breakdown.flight_control    = output_2.wt_flt_ctrl
    output.systems_breakdown.optionals         = output_2.wt_opitems
    output.systems_breakdown.cargo_container   = output_2.wt_cargo_contain
    output.systems_breakdown.crew_bags         = output_2.wt_crew_bags
    output.systems_breakdown.engine_oil        = output_2.wt_eng_oil
    output.systems_breakdown.passenger_service = output_2.wt_pass_serv
    output.systems_breakdown.unusable_fuel     = output_2.wt_unus_fuel
    output.systems_breakdown.air_conditioner   = output_2.wt_ac
    output.systems_breakdown.anti_icing        = output_2.wt_anti_ice
    output.systems_breakdown.wt_fuel_systems   = output_2.wt_fuel_sys
    output.systems_breakdown.wt_misc_systems   = output_2.wt_misc_sys
    
    if propulsors.type == 'electric':
        output.fuel                            = 0
    else:
        output.fuel                            = TOW - output.payload - output.empty

    
    # define weights components
    try: 
        landing_gear_component = vehicle.landing_gear                             # landing gear previously defined
    except AttributeError:                                                      # landing gear not defined
        landing_gear_component = SUAVE.Components.Landing_Gear.Landing_Gear()
        vehicle.landing_gear   = landing_gear_component
    
    control_systems    = SUAVE.Components.Physical_Component()
    electrical_systems = SUAVE.Components.Physical_Component()
    passengers         = SUAVE.Components.Physical_Component()
    furnishings        = SUAVE.Components.Physical_Component()
    air_conditioner    = SUAVE.Components.Physical_Component()
    fuel               = SUAVE.Components.Physical_Component()
    apu                = SUAVE.Components.Physical_Component()
    hydraulics         = SUAVE.Components.Physical_Component()
    instruments        = SUAVE.Components.Physical_Component()
    optionals          = SUAVE.Components.Physical_Component()
    anti_icing         = SUAVE.Components.Physical_Component()
    fuel_systems       = SUAVE.Components.Physical_Component()
    avionics           = SUAVE.Components.Energy.Peripherals.Avionics()
    battery            = SUAVE.Components.Physical_Component()
    battery_payload    = SUAVE.Components.Physical_Component()
    cables             = SUAVE.Components.Physical_Component()
    miscellaneous      = SUAVE.Components.Physical_Component()
    paint              = SUAVE.Components.Physical_Component()
    nacelles           = SUAVE.Components.Physical_Component()
    
    # assign output weights to objects
    landing_gear_component.mass_properties.mass                      = output.landing_gear
    control_systems.mass_properties.mass                             = output.systems_breakdown.flight_control
    electrical_systems.mass_properties.mass                          = output.systems_breakdown.electrical
    passengers.mass_properties.mass                                  = output.payload
    furnishings.mass_properties.mass                                 = output.systems_breakdown.furnish
    air_conditioner.mass_properties.mass                             = output.systems_breakdown.air_conditioner
    fuel.mass_properties.mass                                        = output.fuel
    apu.mass_properties.mass                                         = output.systems_breakdown.apu
    hydraulics.mass_properties.mass                                  = output.systems_breakdown.hydraulics
    avionics.mass_properties.mass                                    = output.systems_breakdown.avionics \
                                                                        + output.systems_breakdown.instruments                  
    optionals.mass_properties.mass                                   = output.systems_breakdown.optionals
    anti_icing.mass_properties.mass                                  = output.systems_breakdown.anti_icing 
    fuel_systems.mass_properties.mass                                = output.systems_breakdown.wt_fuel_systems
    battery.mass_properties.mass                                     = output.battery
    battery_payload.mass_properties.mass                             = output.battery_payload
    cables.mass_properties.mass                                      = output.cables
    miscellaneous.mass_properties.mass                               = output.systems_breakdown.wt_misc_systems
    paint.mass_properties.mass                                       = output.paint
        
    #assign components to vehicle
    vehicle.control_systems                     = control_systems
    vehicle.electrical_systems                  = electrical_systems
    vehicle.avionics                            = avionics
    vehicle.furnishings                         = furnishings
    vehicle.passenger_weights                   = passengers 
    vehicle.air_conditioner                     = air_conditioner
    vehicle.fuel                                = fuel
    vehicle.fuel_systems                        = fuel_systems
    vehicle.apu                                 = apu
    vehicle.anti_icing                          = anti_icing
    vehicle.hydraulics                          = hydraulics
    vehicle.optionals                           = optionals
    vehicle.landing_gear                        = landing_gear_component
    vehicle.battery                             = battery
    vehicle.battery_payload                     = battery_payload
    vehicle.cables                              = cables
    vehicle.paint_w                             = paint
    vehicle.nacelles                            = nacelles


    return output
