## @ingroup Methods-Weights-FLOPS
# systems.py
# 
# Created:  Des 2019, S. Karpuk
# Modified:    

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Units, Data
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Air_conditioning        import  weight_air_cond          as weight_air_cond
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Anti_Icing              import  weight_anti_icing        as weight_anti_icing
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.APU                     import  weight_APU               as weight_APU
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Avionics                import  weight_avionics          as weight_avionics
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Electricals             import  weight_electrical        as weight_electrical
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Furnishings_Equipment   import  weight_furnishings       as weight_firnishings
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Hydraulics              import  weight_hydraulics        as weight_hydraulics
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.Instruments             import  weight_instruments       as weight_instruments
from SUAVE.Methods.Weights.FLOPS.System_and_equipment_items.surface_controls        import  weight_surface_controls  as weight_surface_controls
from SUAVE.Methods.Weights.FLOPS.Operation_items.cargo_containers                   import  weight_cargo_containers  as weight_cargo_containers
from SUAVE.Methods.Weights.FLOPS.Operation_items.crew_baggage                       import  weight_crew_baggage      as weight_crew_baggage
from SUAVE.Methods.Weights.FLOPS.Operation_items.engine_oil                         import  weight_engine_oil        as weight_engine_oil
from SUAVE.Methods.Weights.FLOPS.Operation_items.passenger_service                  import  weight_pass_service      as weight_pass_service
from SUAVE.Methods.Weights.FLOPS.Operation_items.unusable_fuel                      import  weight_unusable_fuel     as weight_unusable_fuel
from SUAVE.Methods.Weights.FLOPS.Propulsion_systems                                 import  weight_fuel_system       as weight_fuel_system

# ----------------------------------------------------------------------
#   Systems
# ----------------------------------------------------------------------

## @ingroup Methods-Weights-FLOPS
def systems(vehicle):
    """ Calculate the weight of different aircraft systems using the FLOPS method 
    
    Assumptions:
        numbers based on FAA regulations and correlations from previous aircraft 

    Source: 
        N/A
                
   Inputs:
       num_seats - total number of seats on the aircraft                                                   [dimensionless]
       ctrl_type - specifies if the control system is fully power, partially powered, or not powered       [dimensionless]
       S_h       - area of the horizontal tail                                                                   [meters**2]
       S_v       - area of the vertical tail                                                                     [meters**2]
       S_gross_w - area of the wing                                                                        [meters**2]
       ac_type   - determines type of instruments, electronics, and operating items based on type of vehicle [dimensionless]
   
   Outputs:
       output             - a data dictionary with fields:
           wt_flt_ctrl    - weight of the flight control system                                               [kilograms]
           wt_apu         - weight of the apu                                                                      [kilograms]
           wt_hyd_pnu     - weight of the hydraulics and pneumatics                                            [kilograms]
           wt_instruments - weight of the instruments and navigational equipment                           [kilograms]
           wt_avionics    - weight of the avionics                                                            [kilograms]
           wt_opitems     - weight of the optional items based on the type of aircraft                         [kilograms]
           wt_elec        - weight of the electrical items                                                        [kilograms]
           wt_ac          - weight of the air conditioning and anti-ice system                                      [kilograms]
           wt_furnish     - weight of the furnishings in the fuselage                                          [kilograms]
       
    Properties Used:
        N/A
    """
    
    # unpack inputs
    vehicle_type = vehicle.type
    Mmax         = vehicle.max_cruise_mach
    range_cr     = vehicle.range
    p_h          = vehicle.cruise_pressure 
    p_0          = vehicle.SL_pressure     
    Nult         = vehicle.envelope.ultimate_load

    TOW = vehicle.mass_properties.max_takeoff
    MZF = vehicle.mass_properties.max_zero_fuel
    
    num_crew   = vehicle.crew
    num_pax    = vehicle.passengers
    num_first  = vehicle.passengers_first
    num_biz    = vehicle.passengers_business
    num_econ   = vehicle.passengers_economy
    wt_cargo   = vehicle.mass_properties.cargo

    fuel_tanks = vehicle.fuel_tanks

    propulsor_name = list(vehicle.propulsors.keys())[0]                     #obtain the key for the propulsor for assignment purposes
    propulsors     = vehicle.propulsors[propulsor_name]
    number_engines = propulsors.number_of_engines
    eng_w          = propulsors.wing_mounted_engines
    eng_f          = propulsors.fuselage_mounted_engines
    engine_d       = propulsors.nacelle_diameter
    APU_tag        = propulsors.APU

    if propulsors.type == 'electric' or propulsors.type == 'cryo-electric' or propulsors.type == 'turboprop':
        eng_thrust = 0
    elif propulsors.type == 'hybrid':
        eng_thrust = propulsors.sealevel_static_thrust_jet
    else:
        eng_thrust = propulsors.sealevel_static_thrust

    S_gross_w    = vehicle.reference_area
    b_w          = vehicle.wings['main_wing'].spans.projected
    sweep_w      = vehicle.wings['main_wing'].sweeps.quarter_chord
    NBAYS        = vehicle.wings['main_wing'].cabin.number_of_bays  
    VARSWP       = vehicle.wings['main_wing'].sweeps.variable
    sweep_f      = vehicle.wings['main_wing'].cabin_sweep_leading_edge
    
    FPAREA       = vehicle.fuselages['fuselage'].areas.top_projected
    cabin_area   = vehicle.fuselages['fuselage'].areas.cabin_floor
    CARGF        = vehicle.fuselages['fuselage'].CARGF
    n_fus        = vehicle.fuselages['fuselage'].number_of_fuselages
    l_fus        = vehicle.fuselages['fuselage'].lengths.total
    l_cab        = vehicle.fuselages['fuselage'].lengths.cabin 
    h_fus        = vehicle.fuselages['fuselage'].heights.maximum
    w_fus        = vehicle.fuselages['fuselage'].width

    w_misc       = vehicle.miscellaneous.weight
    
    
    # process   
    # APU Group ----------------------------------------------------------------------------------------------  
    if APU_tag == True:
            apu_wt = weight_APU(num_pax,FPAREA,n_fus)
    else:
        apu_wt = 0.0 * Units.lb     #   no apu if less than 9 seats
    #---------------------------------------------------------------------------------------------------------

    # Avionics Group -----------------------------------------------------------------------------------------
    avionics_wt = weight_avionics(vehicle_type,range_cr,num_crew,FPAREA,n_fus,l_fus,h_fus,Mmax,0)
    #---------------------------------------------------------------------------------------------------------

    # Electrical Group ---------------------------------------------------------------------------------------
    electrical_wt = weight_electrical(vehicle_type,l_fus,w_fus,n_fus,number_engines, \
                                      num_crew,num_pax,Mmax,b_w)
    #---------------------------------------------------------------------------------------------------------

    # Furnishings Group --------------------------------------------------------------------------------------   
    furnish_wt = weight_firnishings(vehicle_type,num_crew,num_first,num_biz,num_econ,l_cab,n_fus,Mmax, \
                                    l_fus,CARGF,cabin_area,w_fus,h_fus,NBAYS,sweep_f)
    #---------------------------------------------------------------------------------------------------------
    
    # Hydraulics & Pneumatics Group --------------------------------------------------------------------------
    hyd_pnu_wt = weight_hydraulics(FPAREA,Mmax,S_gross_w,eng_w,eng_f,VARSWP,n_fus)
    #---------------------------------------------------------------------------------------------------------

    # Instruments Group --------------------------------------------------------------------------------------
    instruments_wt = weight_instruments(vehicle_type,FPAREA,Mmax,num_crew,eng_w,eng_f,l_fus,h_fus,n_fus)
    #---------------------------------------------------------------------------------------------------------

    # Surface_controls Group ---------------------------------------------------------------------------------
    SFLAP        = vehicle.wings['main_wing'].areas.flap 
    flt_ctrl_wt  = weight_surface_controls(vehicle_type,TOW,Mmax,SFLAP,p_0,p_h,Nult,S_gross_w)
    #---------------------------------------------------------------------------------------------------------

    # Operating items Group ----------------------------------------------------------------------------------
    # Cargo containers
    if vehicle.fuselages['fuselage'].cargo_container == True:
        cargo_contain_wt = weight_cargo_containers(wt_cargo)
    else:
        cargo_contain_wt = 0.0


    # Crew baggage
    crew_bags_wt = weight_crew_baggage(num_pax,vehicle_type)

    if propulsors.type != 'electric' or propulsors.type != 'propeller electric' or propulsors.type != 'cryo-electric' or propulsors.type != 'propeller cryo-electric':
        if propulsors.type == 'hybrid':
            
            if propulsors.thrust.hybridization != 0:
                number_engines  = propulsors.jet_engines
                number_wing_eng = propulsors.wing_mounted_jet_engines
                number_fuse_eng = propulsors.fuselage_mounted_jet_engines

                # Engine_oil
                eng_oil_wt = weight_engine_oil(eng_thrust,number_engines)
    
                # Unusable_fuel
                unus_fuel_wt = weight_unusable_fuel(vehicle_type,number_engines,eng_thrust,S_gross_w, \
                                                fuel_tanks,TOW-MZF)

                # Hydraulics & Pneumatics Group 
                hyd_pnu_wt = weight_hydraulics(FPAREA,Mmax,S_gross_w,number_wing_eng,number_fuse_eng,VARSWP)

                # Fuel system Group 
                fuel_sys_wt = weight_fuel_system(number_wing_eng,number_fuse_eng,vehicle_type,TOW-MZF,Mmax,fuel_tanks)

                # Electrical Group 
                electrical_wt = weight_electrical(vehicle_type,l_fus,w_fus,n_fus,number_engines, \
                                                    num_crew,num_pax,Mmax,b_w)

                # Instruments Group 
                instruments_wt = weight_instruments(vehicle_type,FPAREA,Mmax,num_crew,number_wing_eng,number_fuse_eng,l_fus,h_fus,n_fus)

            else:
                # Engine_oil
                eng_oil_wt = weight_engine_oil(eng_thrust,number_engines)
    
                # Unusable_fuel
                unus_fuel_wt = weight_unusable_fuel(vehicle_type,number_engines,eng_thrust,S_gross_w, \
                                                        fuel_tanks,TOW-MZF)

                # Hydraulics & Pneumatics Group 
                hyd_pnu_wt = weight_hydraulics(FPAREA,Mmax,S_gross_w,eng_w,eng_f,VARSWP)

                # Fuel system Group 
                fuel_sys_wt = weight_fuel_system(eng_w,eng_f,vehicle_type,TOW-MZF,Mmax,fuel_tanks)

                # Electrical Group 
                electrical_wt = weight_electrical(vehicle_type,l_fus,w_fus,n_fus,number_engines, \
                                                    num_crew,num_pax,Mmax,b_w)

                # Instruments Group 
                instruments_wt = weight_instruments(vehicle_type,FPAREA,Mmax,num_crew,eng_w,eng_f,l_fus,h_fus,n_fus)
        else:
            # Engine_oil
            eng_oil_wt = weight_engine_oil(eng_thrust,number_engines)
    
            # Unusable_fuel
            unus_fuel_wt = weight_unusable_fuel(vehicle_type,number_engines,eng_thrust,S_gross_w, \
                                                fuel_tanks,TOW-MZF)

            # Hydraulics & Pneumatics Group 
            hyd_pnu_wt = weight_hydraulics(FPAREA,Mmax,S_gross_w,eng_w,eng_f,VARSWP)

            # Fuel system Group 
            fuel_sys_wt = weight_fuel_system(eng_w,eng_f,vehicle_type,TOW-MZF,Mmax,fuel_tanks)

            # Electrical Group 
            electrical_wt = weight_electrical(vehicle_type,l_fus,w_fus,n_fus,number_engines, \
                                      num_crew,num_pax,Mmax,b_w)
                    
            # Instruments Group 
            instruments_wt = weight_instruments(vehicle_type,FPAREA,Mmax,num_crew,eng_w,eng_f,l_fus,h_fus,n_fus)

    else:
        # Engine_oil
        eng_oil_wt   = 0
        
        # Unusable_fuel
        unus_fuel_wt = 0

        # Hydraulics & Pneumatics Group 
        hyd_pnu_wt   = 0

        # Fuel system Group 
        fuel_sys_wt  = 0

    # Passenger_service
    pass_serv_wt = weight_pass_service(num_first,num_biz,num_econ,range_cr,Mmax)


    opitems_wt = cargo_contain_wt + crew_bags_wt + eng_oil_wt + pass_serv_wt + unus_fuel_wt
    #---------------------------------------------------------------------------------------------------------
       
    # Environmental Control ----------------------------------------------------------------------------------
    ac_indicator = vehicle.fuselages['fuselage'].air_conditioning
    if ac_indicator != True:
        ac_wt = 0
    else:           
        ac_wt = weight_air_cond(vehicle_type,FPAREA,h_fus,vehicle.passengers,Mmax,avionics_wt,number_engines, \
                                eng_thrust/number_engines,n_fus)
    #---------------------------------------------------------------------------------------------------------

    # Anti-icing ---------------------------------------------------------------------------------------------
    if vehicle_type != 'transport':
        anti_ice_wt = 0
    else:
        anti_ice_wt = weight_anti_icing(b_w,sweep_w,engine_d,number_engines,w_fus)
    #---------------------------------------------------------------------------------------------------------

    # Fuel system Group --------------------------------------------------------------------------------------
    fuel_sys_wt = weight_fuel_system(eng_w,eng_f,vehicle_type,TOW-MZF,Mmax,fuel_tanks)
    #---------------------------------------------------------------------------------------------------------


    # Miscelaneous system Group ------------------------------------------------------------------------------
    misc_sys_wt = w_misc
    #---------------------------------------------------------------------------------------------------------

    # packup outputs
    output = Data()
    output.wt_apu           = apu_wt
    output.wt_avionics      = avionics_wt
    output.wt_electrical    = electrical_wt
    output.wt_furnish       = furnish_wt
    output.wt_hyd_pnu       = hyd_pnu_wt
    output.wt_instruments   = instruments_wt
    output.wt_flt_ctrl      = flt_ctrl_wt
    output.wt_opitems       = opitems_wt
    output.wt_cargo_contain = cargo_contain_wt
    output.wt_crew_bags     = crew_bags_wt
    output.wt_eng_oil       = eng_oil_wt
    output.wt_pass_serv     = pass_serv_wt
    output.wt_unus_fuel     = unus_fuel_wt
    output.wt_ac            = ac_wt
    output.wt_anti_ice      = anti_ice_wt
    output.wt_fuel_sys      = fuel_sys_wt
    output.wt_misc_sys      = misc_sys_wt

    output.wt_systems       = output.wt_apu + output.wt_avionics + output.wt_electrical + output.wt_furnish + \
                              output.wt_hyd_pnu + output.wt_instruments + output.wt_flt_ctrl + output.wt_opitems  + \
                              output.wt_ac + output.wt_anti_ice + output.wt_fuel_sys + output.wt_misc_sys

    
    return output
