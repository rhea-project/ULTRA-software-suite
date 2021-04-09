## @ingroup Methods-Costs-Operating_Costs
# compute_operating_costs.py
#
# Created:  

# suave imports
import numpy as np
from SUAVE.Core import Units,Data

# ----------------------------------------------------------------------
#  Compute operating costs
# ----------------------------------------------------------------------

## @ingroup Methods-Costs-Operating_Costs
def compute_operating_costs(vehicle):
    """ This function calculated Direct Operating Costs (DOC) using the AEA 1989 method

        Assumptions:
            take-off propulsive efficiency of 65% for propeller- driven aircraft
    
        Source:

            J. Hoelzen, Y. Liu, B. Bensmann, C. Winnefeld, A. Elham, J. Friedrichs, and R. Hanke-Rauschenbach
                'Conceptual Design of Operation Strategies for Hybrid Electric Aircraft'
            S. Gudmundsson, 'General Aviation Aircraft Design. Applied Methods and Procedures'


        Inputs:
            vehicle   
    
        Outputs:
            DOC
    
        Properties Used:
            N/A
        """    
    
    # Unpack inputs
    #-----------------------------------------------------------------------------
    costs   = vehicle.costs.operating

    CPI     = vehicle.costs.CPI 
    FC      = costs.flight_cycles
    LR      = costs.labor_costs              # Labor costs in EUR/hr
    p_fuel  = costs.fuel_price 
    k_Rep   = costs.repair_costs             # Repair costs per flight
    kNav    = costs.navigation_fees
    kLand   = costs.landing_fees
    kGround = costs.landing_ground   
    CC      = costs.crew_complements
    sFC     = costs.pilot_rate
    sFA     = costs.crew_rate
    nFC     = vehicle.pilots
    nFA     = vehicle.flight_attendants

    Range   = vehicle.range / 1000
    t_Total = vehicle.mission_time

    DP      = costs.depreciate_years
    IR      = costs.interest_rate 
    fins    = costs.insure_rate 

    MTOW    = vehicle.mass_properties.max_takeoff 
    OEW     = vehicle.mass_properties.operating_empty
    Wp      = vehicle.mass_properties.payload
    Wengine = vehicle.mass_properties.engine
    FB      = vehicle.mass_properties.fuel_burn

    b = vehicle.wings['main_wing'].spans.projected
    l = vehicle.fuselages['fuselage'].lengths.total 


    # Calcualte costs for a General Aviation aircraft 
    #-----------------------------------------------------------------------------
    if vehicle.type == 'general aviation' or vehicle.type == 'general aviation business':

        # Unpack inputs
        #---------------------------------------------------
        R_stor = costs.storage_rate                                         # Storage rate per month
        R_insp = costs.inspection
        C_AC   = costs.aircraft_insured_value 
        pAC    = costs.capital.price.aircraft

        Neng   = vehicle.propulsors.turbofan_hybrid.number_of_engines


        # Constants
        #---------------------------------------------------
        # Note: Refer to Gudmundsson for a detailed description for the F-matrix (Eq 2-28) 
        # and adjust it according to your aircraft 
        if vehicle.type == 'general aviation':
            F = [0.3, 0, 0.02, 0.02, 0.02, 0.04, 0.01, 0.02, 0]
        else:
            F = [2.0, 0.2, 0.2, 0.2, 0.1, 0.2, 0.5]

        # Maintenance Costs
        #---------------------------------------------------
        C_AP = np.sum(F) * LR * t_Total * FC

        # Storage Costs
        #---------------------------------------------------
        C_STOR = 12 * R_stor

        # Fuel Costs
        #---------------------------------------------------
        C_FUEL = FB * t_Total * FC * p_fuel

        # Insurance Costs
        #---------------------------------------------------
        C_INS = 500 + 0.015 * C_AC

        # Annual Inspection Costs
        #---------------------------------------------------
        C_INSP = R_insp

        # Engine Overhaul
        #---------------------------------------------------
        if vehicle.type == 'general aviation':
            C_OVER = 5 * Neng * t_Total * FC 
        else:
            C_OVER = 7.5 * Neng * t_Total * FC 

        # Crew costs
        #---------------------------------------------------
        if vehicle.type != 'general aviation':
            C_Crew = CC * (sFA*nFA + sFC*nFC)
        else:
            C_Crew = 0.0

        # Loan Costs
        #---------------------------------------------------       
        C_LOAN = 12 * pAC / (1-1/(1+IR)**(12*DP))

        DOCtot =  (C_AP + C_STOR + C_FUEL + C_INS + C_INSP + C_OVER + C_LOAN + C_Crew)/FC     # DOC per flight

        # Pack outputs
        DOC = Data()
        DOC.Fuel_Energy    = C_FUEL / FC
        DOC.Crew           = C_Crew/FC
        DOC.Maintenance    = (C_AP+C_INSP)/FC
        DOC.Capital        = (C_LOAN+C_INS)/FC
        DOC.Total          = DOCtot
        DOC.flight_cycles  = FC


    # Calculate costs for a battery hybrid aircraft
    #-----------------------------------------------------------------------------

    else:
        # Unpack unputs
        #---------------------------------------------------
        fRV     = costs.capital.residual_value

        k_Ma    = costs.maintenance_gain

        kS_AF   = costs.capital.spare_parts.airframe
        kS_GT   = costs.capital.spare_parts.gas_turbine

        pAF     = costs.capital.price.airframe
        
        Wbat    = 0
 

        if vehicle.propulsion_type == 'battery-hybrid' and vehicle.type == 'transport':
            p_elec  = costs.electricity_price

            fRV_Bat = costs.capital.residual_value_battery
            pEM     = costs.capital.price.electric_motor
            pPMAD   = costs.capital.price.PMAD
            pBat    = costs.capital.price.Battery

            N_Bat_c = costs.capital.battery_cycles

            kS_EM   = costs.capital.spare_parts.electric_motor
            kS_PMAD = costs.capital.spare_parts.kS_PMAD

            Hp       = vehicle.propulsors.turbofan_hybrid.thrust.hybridization
            Neng     = vehicle.propulsors.turbofan_hybrid.number_of_engines
            T_SL     = vehicle.propulsors.turbofan_hybrid.sealevel_static_thrust_jet / (1 - Hp) 
            P_EM     = vehicle.propulsors.turbofan_hybrid.motor.power/1000
            eta_EM   = vehicle.propulsors.turbofan_hybrid.motor.motor_efficiency 
            eta_PMAD = vehicle.propulsors.turbofan_hybrid.inverter.efficiency

            Ebat    = vehicle.battery_energy*0.00000027778      # kWh
            Wbat    = vehicle.mass_properties.battery

            # Capital costs of batteries
            DP_Bat      = 3 * N_Bat_c/FC
            aBat        = IR * (1-fRV_Bat*(1/(1+IR))**DP_Bat)/(1-(1/(1+IR))**DP_Bat)
            DOC_CAP_Bat = 3 * Ebat * pBat * (aBat + fins) * CPI

        elif (vehicle.propulsion_type == 'cryo-electric' or vehicle.propulsion_type == 'all-electric') and vehicle.type == 'transport' :

            Neng     = vehicle.propulsors.network.number_of_engines

            p_elec  = costs.electricity_price

            fRV_Bat = costs.capital.residual_value_battery
            pEM     = costs.capital.price.electric_motor
            pPMAD   = costs.capital.price.PMAD
            pBat    = costs.capital.price.Battery

            Ebat    = vehicle.battery_energy*0.00000027778      # kWh
            Wbat    = vehicle.mass_properties.battery

            N_Bat_c = costs.capital.battery_cycles

            kS_EM   = costs.capital.spare_parts.electric_motor
            kS_PMAD = costs.capital.spare_parts.kS_PMAD

            P_EM     = vehicle.propulsors.network.motor.sea_level_power * Neng
            eta_EM   = vehicle.propulsors.network.motor.motor_efficiency 
            eta_PMAD = vehicle.propulsors.network.inverter.efficiency

            V1 = vehicle.liftoff_speed
            Hp = 1.0

            # Capital costs of batteries
            DP_Bat      = 3 * N_Bat_c/FC
            aBat        = IR * (1-fRV_Bat*(1/(1+IR))**DP_Bat)/(1-(1/(1+IR))**DP_Bat)
            DOC_CAP_Bat = 3 * Ebat * pBat * (aBat + fins) * CPI


        elif vehicle.propulsion_type == 'turboprop' and vehicle.type == 'transport' :

            Neng     = vehicle.propulsors.turboprop.number_of_engines
            P_GT     = vehicle.propulsors.turboprop.engine.sea_level_power/1000 * Neng
            pGT      = costs.capital.price.gas_turbine                               # EUR/kW
            
            V1 = vehicle.liftoff_speed

            Hp          = 0.0
            pEM         = 0.0
            N_Bat_c     = 0.0
            pPMAD       = 0.0
            pBat        = 0.0
            P_EM        = 0.0
            fRV_Bat     = 0.0
            kS_EM       = 0.0
            kS_PMAD     = 0.0
            eta_EM      = 0.5
            eta_PMAD    = 0.5
            Ebat        = 0.0
            Wbat        = 0.0
            DOC_CAP_Bat = 0.0
            p_elec      = 0.0

        else:
            Neng        = vehicle.propulsors.turbofan.number_of_engines
            T_SL        = vehicle.propulsors.turbofan.sealevel_static_thrust
            Hp          = 0.0
            pEM         = 0.0
            N_Bat_c     = 0.0
            pPMAD       = 0.0
            pBat        = 0.0
            P_EM        = 0.0
            fRV_Bat     = 0.0
            kS_EM       = 0.0
            kS_PMAD     = 0.0
            eta_EM      = 0.5
            eta_PMAD    = 0.5
            Ebat        = 0.0
            Wbat        = 0.0
            DOC_CAP_Bat = 0.0
            p_elec      = 0.0

        # Energy costs
        #---------------------------------------------------
        if (vehicle.propulsion_type == 'battery-hybrid' or vehicle.propulsion_type == 'cryo-electric' or vehicle.propulsion_type == 'cryo-electric') and vehicle.type == 'transport' :
            DOC_Energy = FC * (FB*p_fuel + Ebat*p_elec)
        else:
            DOC_Energy = FC * FB*p_fuel

        # Crew costs
        #---------------------------------------------------
        DOC_Crew   = CC * (sFA*nFA + sFC*nFC) * CPI

        # Maintenance costs
        #---------------------------------------------------

        # Airframe material costs
        Waf        = OEW - Wengine - Wbat
        DOC_AF_mat = (Waf/1000 * (0.21*t_Total + 13.7)+57.5 ) * CPI

        # Airframe personnel costs
        B = 2
        DOC_AF_per = LR * ((0.01*Waf/1000 + 0.655)*t_Total + (0.01*Waf/1000 + 0.254)) * CPI * (1+B)

        # Engine maintenance costs
        if Hp == 1.0:
            eta = 0.91
        else:
            eta = 1.0

        if (vehicle.propulsion_type == 'battery-hybrid' or vehicle.propulsion_type == 'cryo-electric' or vehicle.propulsion_type == 'cryo-electric') and vehicle.type == 'transport' :    
            DOC_Eng  = eta * (7.621e-4*0.65*P_EM/V1 + 30.5*t_Total + 10.6) * Neng * CPI

        elif vehicle.propulsion_type == 'turboprop':
            DOC_Eng  = eta * (7.621e-4*0.65*P_GT/V1 + 30.5*t_Total + 10.6) * Neng * CPI

        else:
            DOC_Eng  = eta * (1.5*T_SL/9.81/1000 + 30.5*t_Total + 10.6) * Neng * CPI

        # Technology costs
        DOC_Tech = 0.0 #5000 * (b*l)**0.75 

        DOC_Ma = (DOC_AF_mat+DOC_AF_per+DOC_Eng+DOC_Tech) * FC * k_Ma * CPI

        # Depreciation costs
        #---------------------------------------------------

        # Capital costs without batteries
        a          = IR * (1-fRV*(1/(1+IR))**DP)/(1-(1/(1+IR))**DP)

        if vehicle.propulsion_type == 'turboprop':
            PGT = P_GT * pGT
        elif  vehicle.propulsion_type == 'cryo-electric' or vehicle.propulsion_type == 'all-electric':
            PGT = 0.0
        else:
            PGT = 2500 * Wengine        # for jet engines only

        DOC_Cap_AC = ((pAF*Waf*(1+kS_AF)+PGT*(1+kS_GT)+pEM*P_EM*(1+kS_EM)+ \
                      pPMAD*(P_EM/eta_EM/eta_PMAD)*(1+kS_PMAD))*(a+fins))


        DOC_Cap = (DOC_Cap_AC + DOC_CAP_Bat) * CPI 

        # Fees costs
        #---------------------------------------------------

        # Landing
        DOC_LG = kLand*MTOW*FC * CPI 

        # Ground operations
        DOC_Ground = (kGround*Wp)*FC * CPI 

        # Navigation
        DOC_Nav = (kNav * Range/100 * np.sqrt(MTOW/50000))*FC * CPI 

        DOC_Fees = DOC_LG + DOC_Ground + DOC_Nav

        DOCtot =  (DOC_Energy + DOC_Crew + DOC_Ma + DOC_Cap + DOC_Fees)/FC     # DOC per flight   

        # Pack outputs
        DOC = Data()
        DOC.Energy         = DOC_Energy/FC
        DOC.Fuel_Energy    = FB*p_fuel
        DOC.Battery_Energy = Ebat*p_elec
        DOC.Crew           = DOC_Crew/FC
        DOC.Maintenance    = DOC_Ma/FC
        DOC.Capital        = DOC_Cap/FC
        DOC.Fees           = DOC_Fees/FC
        DOC.Total          = DOCtot
        DOC.flight_cycles  = FC

    
    return DOC
