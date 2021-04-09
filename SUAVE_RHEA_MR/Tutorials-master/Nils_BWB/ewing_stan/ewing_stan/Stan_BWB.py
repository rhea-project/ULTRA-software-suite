# tut_mission_B737_AVL.py
# 
# Created:  Mar 2018, SUAVE Team

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

# Python Imports
import numpy as np
import pylab as plt

# SUAVE Imports
import SUAVE
from SUAVE.Core import Units, Data
import numpy as np
from SUAVE.Components.Energy.Networks.Solar import Solar
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Methods.Power.Battery.Sizing import initialize_from_energy_and_power, initialize_from_mass


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    configs, analyses = full_setup()

    simple_sizing(configs)

    configs.finalize()
    analyses.finalize()

    # weight analysis
    weights = analyses.configs.base.weights
    breakdown = weights.evaluate()      

    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate()

    # plt the old results
    plot_mission(results)

    return

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():

    # vehicle data
    vehicle  = vehicle_setup()
    configs  = configs_setup(vehicle)

    # vehicle analyses
    configs_analyses = analyses_setup(configs)

    # mission analyses
    mission  = mission_setup(configs_analyses)
    missions_analyses = missions_setup(mission)

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses

    return configs, analyses

# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):

    analyses = SUAVE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis

    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = SUAVE.Analyses.Vehicle()

    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = SUAVE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)

    # ------------------------------------------------------------------
    #  Weights
    weights = SUAVE.Analyses.Weights.Weights_Tube_Wing()
    weights.vehicle = vehicle
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.AVLQ3D()
    #aerodynamics = SUAVE.Analyses.Aerodynamics.AVL()
    aerodynamics.settings.spanwise_vortices     = 25 
    aerodynamics.settings.chordwise_vortices    = 10
    aerodynamics.geometry = vehicle
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = SUAVE.Analyses.Stability.AVL()
    #stability.settings.spanwise_vortex_density                  = 3
    stability.geometry = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= SUAVE.Analyses.Energy.Energy()
    energy.network = vehicle.propulsors 
    analyses.append(energy)

    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = SUAVE.Analyses.Planets.Planet()
    analyses.append(planet)

    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere)   

    return analyses    

# ----------------------------------------------------------------------
#   Define the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = SUAVE.Vehicle()
    vehicle.tag = 'Boeing_737-800'    
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    

    # mass properties
    vehicle.mass_properties.takeoff             = 580. * Units.kg
    #vehicle.mass_properties.operating_empty    = 385. * Units.kg
    vehicle.mass_properties.max_takeoff         = 600. * Units.kg
    vehicle.mass_properties.payload             = 200. * Units.kg  
    
    # envelope properties
    #vehicle.envelope.ultimate_load = 2.5
    #vehicle.envelope.limit_load    = 1.5

    # basic parameters
    vehicle.reference_area         = 16.5 * Units['meters**2']  
    vehicle.passengers             = 170


    # ------------------------------------------------------------------        
    #  Landing Gear
    # ------------------------------------------------------------------        
    # used for noise calculations
    landing_gear = SUAVE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag = "main_landing_gear"
    
    landing_gear.main_tire_diameter = 0.25000 * Units.m
    landing_gear.nose_tire_diameter = 0.20 * Units.m
    landing_gear.main_strut_length  = 0.3 * Units.m
    landing_gear.nose_strut_length  = 0.6 * Units.m
    landing_gear.main_units  = 2    #number of main landing gear units
    landing_gear.nose_units  = 1    #number of nose landing gear
    landing_gear.main_wheels = 1    #number of wheels on the main landing gear
    landing_gear.nose_wheels = 1    #number of wheels on the nose landing gear
    vehicle.landing_gear = landing_gear

    # ------------------------------------------------------------------        
    #   Main Wing
    # ------------------------------------------------------------------        
    
    wing = SUAVE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'

    wing.areas.reference         = vehicle.reference_area
    wing.aspect_ratio            = (wing.spans.projected**2.) / wing.areas.reference
    wing.spans.projected         = 11. * Units.meter
	
    wing.chords.root             = 3.4 * Units.meter
    wing.chords.seat             = 3.2 * Units.meter
    wing.chords.kink             = 1.55 * Units.meter
    wing.chords.tip              = 0.4 * Units.meter
    wing.chords.mid_wing         = (wing.chords.kink + wing.chords.tip) / 2. 

    spanwise_locations_root        = 0.0
    spanwise_locations_seat        = 0.6 * Units.meter / wing.spans.projected * 2
    spanwise_locations_kink        = 1.5 * Units.meter / wing.spans.projected * 2
    spanwise_locations_tip         = 0.97
    spanwise_locations_mid_wing    = (spanwise_locations_kink + spanwise_locations_tip) / 2

    wing.areas.reference         = ((wing.chords.root + wing.chords.seat) * spanwise_locations_seat +
                                    (wing.chords.seat + wing.chords.kink) * (spanwise_locations_kink - spanwise_locations_seat) +
                                    (wing.chords.tip + wing.chords.kink) * (spanwise_locations_tip - spanwise_locations_kink)) * wing.spans.projected / 2.
	

    vehicle.reference_area       = wing.areas.reference

    wing.taper                   = wing.chords.tip / wing.chords.root
    wing.chords.mean_aerodynamic = wing.areas.reference / wing.spans.projected
	
    wing.twists.root             = -2.5 * Units.degrees
    wing.twists.seat             = -2.0 * Units.degrees
    wing.twists.kink             = 0.0 * Units.degrees
    wing.twists.tip              = -7.0 * Units.degrees
    wing.twists.mid_wing         = (wing.twists.tip - wing.twists.kink) / 4.
	
    wing.origin                  = [0.,0.,0.] # meters
    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = True
    wing.dynamic_pressure_ratio  = 1.0
    wing.manual_segments         = False


    segment = SUAVE.Components.Wings.Segment()
    
    segment.tag                   = 'root'
    segment.percent_span_location = spanwise_locations_root
    segment.twist                 = wing.twists.root
    segment.root_chord_percent    = 1
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 8.0 * Units.degrees
	
    section_1_airfoil = SUAVE.Components.Wings.Airfoils.Airfoil()
    section_1_airfoil.coordinate_file = 'NACA65_410.dat'   
    segment.append_airfoil(section_1_airfoil)
    wing.Segments.append(segment) 
	
    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'seat'
    segment.percent_span_location = spanwise_locations_seat
    segment.twist                 = wing.twists.seat
    segment.root_chord_percent    = wing.chords.seat / wing.chords.root
    segment.dihedral_outboard     = 1.   * Units.degrees
    segment.sweeps.quarter_chord  = 7.0 * Units.degrees
    section_2_airfoil = SUAVE.Components.Wings.Airfoils.Airfoil()
    section_2_airfoil.coordinate_file = 'NACA65_410.dat'
    segment.append_airfoil(section_2_airfoil)
    wing.Segments.append(segment)   
    
    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'kink'
    segment.percent_span_location = spanwise_locations_kink
    segment.twist                 = wing.twists.kink
    segment.root_chord_percent    = wing.chords.kink / wing.chords.root
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 25 * Units.degrees
    section_3_airfoil = SUAVE.Components.Wings.Airfoils.Airfoil()
    section_3_airfoil.coordinate_file = 'NACA65_410.dat'
    segment.append_airfoil(section_3_airfoil)
    wing.Segments.append(segment)

    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'tip'
    segment.percent_span_location = spanwise_locations_tip
    segment.twist                 = wing.twists.tip
    segment.root_chord_percent    = wing.chords.tip / wing.chords.root
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 25. * Units.degrees
    section_4_airfoil = SUAVE.Components.Wings.Airfoils.Airfoil()
    section_4_airfoil.coordinate_file = 'NACA65_410.dat'
    segment.append_airfoil(section_4_airfoil)
    #segment.thickness_to_chord    = 0.175
    wing.Segments.append(segment)

    
    # ------------------------------------------------------------------
    #   Flaps
    # ------------------------------------------------------------------
    wing.flaps.chord      =  0.30   # 30% of the chord
    wing.flaps.span_start =  0.10   # 10% of the span
    wing.flaps.span_end   =  0.75
    wing.flaps.type       = 'double_slotted'

    # add to vehicle
    vehicle.append_component(wing)

    
    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    
    wing = SUAVE.Components.Wings.Wing()
    wing.tag = 'vertical_stabilizer'    

    wing.sweeps.quarter_chord    = 20. * Units.deg
    wing.thickness_to_chord      = 0.12
    wing.spans.projected         = 0.5
    wing.chords.root             = 1.5 * Units.meter
    wing.chords.tip              = 1.0 * Units.meter
    wing.taper                   = wing.chords.tip / wing.chords.root
    wing.areas.reference         = (wing.chords.root + wing.chords.tip) / 2 * wing.spans.projected
    wing.aspect_ratio            = wing.spans.projected * wing.spans.projected / wing.areas.reference
    wing.chords.mean_aerodynamic = 1.0 * Units.meter
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [2.2, 0.0, 0.0]  # meters
    wing.vertical                = True 
    wing.symmetric               = False
    wing.t_tail                  = False
    wing.dynamic_pressure_ratio  = 1.0
        
    # add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    """
    fuselage = SUAVE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'
    
    fuselage.number_coach_seats    = vehicle.passengers
    fuselage.seats_abreast         = 6
    fuselage.seat_pitch            = 1     * Units.meter
    fuselage.fineness.nose         = 1.6
    fuselage.fineness.tail         = 2.
    fuselage.lengths.nose          = 6.4   * Units.meter
    fuselage.lengths.tail          = 8.0   * Units.meter
    fuselage.lengths.cabin         = 28.85 * Units.meter
    fuselage.lengths.total         = 38.02 * Units.meter
    fuselage.lengths.fore_space    = 6.    * Units.meter
    fuselage.lengths.aft_space     = 5.    * Units.meter
    fuselage.width                 = 3.74  * Units.meter
    fuselage.heights.maximum       = 3.74  * Units.meter
    fuselage.effective_diameter    = 3.74     * Units.meter
    fuselage.areas.side_projected  = 142.1948 * Units['meters**2'] 
    fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
    fuselage.areas.front_projected = 12.57    * Units['meters**2'] 
    fuselage.differential_pressure = 5.0e4 * Units.pascal # Maximum differential pressure
    
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter
    
    # add to vehicle
    vehicle.append_component(fuselage)"""

   #------------------------------------------------------------------
    # Propulsor
    #------------------------------------------------------------------
    
    # build network
    net = Solar()
    net.number_of_engines = 1.
    net.nacelle_diameter  = 0.2 * Units.meters
    net.engine_length     = 0.1 * Units.meters
    net.areas             = Data()
    net.areas.wetted      = 0.1*(2.*np.pi*0.2/2.)
    
    # Component 1 the Sun?
    sun = SUAVE.Components.Energy.Processes.Solar_Radiation()
    net.solar_flux = sun
    
    # Component 2 the solar panels
    panel = SUAVE.Components.Energy.Converters.Solar_Panel()
    panel.area                 = (vehicle.reference_area - 3.0) * 0.9
    panel.efficiency           = 0.224
    panel.mass_properties.mass = panel.area*(0.60 * Units.kg)
    net.solar_panel            = panel
    
    # Component 3 the ESC
    # Emsiso EmDrive 500 Controller
    esc = SUAVE.Components.Energy.Distributors.Electronic_Speed_Controller()
    esc.efficiency = 0.97 # Gundlach for brushless motors
    esc.mass_properties.mass = 4.9  * Units.kg
    esc.mass_properties.center_of_gravity         = np.array([1.4,0.0,-0.4])
    
    net.esc        = esc

    # Component 5 the Propeller
    # Design the Propeller
    prop_attributes = Data()
    prop_attributes.number_blades       = 2.0
    prop_attributes.freestream_velocity = 54.0 * Units['m/s'] # freestream
    prop_attributes.angular_velocity    = 1250. * Units['rpm']
    prop_attributes.tip_radius          = 0.875 * Units.meters
    prop_attributes.hub_radius          = 0.08 * Units.meters
    prop_attributes.design_Cl           = 0.75
    prop_attributes.design_altitude     = 2.0 * Units.km
    prop_attributes.design_thrust       = 0.
    prop_attributes.design_power        = 14000.0 * Units.watts
    prop_attributes                     = propeller_design(prop_attributes)
    
    prop = SUAVE.Components.Energy.Converters.Propeller()
    prop.prop_attributes = prop_attributes


    '''prop.materials = SUAVE.Core.Container()
    prop.materials.skin_materials = SUAVE.Core.Container()
    prop.materials.spar_materials = SUAVE.Core.Container()
    prop.materials.flap_materials = SUAVE.Core.Container()
    prop.materials.skin_materials = SUAVE.Core.Container()
    prop.materials.rib_materials = SUAVE.Core.Container()
    prop.materials.root_materials = SUAVE.Core.Container()

    prop.materials.skin_materials.torsion_carrier       = bid_carbon
    prop.materials.spar_materials.shear_carrier         = bid_carbon
    prop.materials.flap_materials.bending_carrier       = ud_carbon
    prop.materials.skin_materials.core                  = core
    prop.materials.rib_materials.structural             = aluminium_rib
    # prop.materials.ribMat.minimum_width       
    prop.materials.root_materials.structural            = ud_carbon
    prop.materials.skin_materials.leading_edge          = ud_carbon
    prop.materials.skin_materials.adhesive              = epoxy
    prop.materials.skin_materials.cover                 = epoxy'''


    prop.mass_properties.mass = 3. * Units.kg
    prop.mass_properties.center_of_gravity         = np.array([3.55,0.0,0.31])

    net.propeller        = prop

    # Component 4 the Motor
    # values for Emrax 228 MV motor
    motor = SUAVE.Components.Energy.Converters.Motor()
    motor.resistance           = 18.0e-3
    motor.no_load_current      = 4.17 * Units.ampere
    motor.speed_constant       = 11.0 * Units['rpm'] # RPM/volt converted to (rad/s)/volt    
    motor.propeller_radius     = prop.prop_attributes.tip_radius
    motor.propeller_Cp         = prop.prop_attributes.Cp
    motor.gear_ratio           = 1. # Gear ratio
    motor.gearbox_efficiency   = 1. # Gear box efficiency
    motor.expected_current     = 120. # Expected current
    motor.mass_properties.mass = 12.3  * Units.kg
    motor.mass_properties.center_of_gravity         = np.array([3.35,0.0,0.31])
    net.motor                  = motor    
    
    # Component 6 the Payload
    payload = SUAVE.Components.Energy.Peripherals.Payload()
    payload.power_draw           = 0. * Units.watts 
    payload.mass_properties.mass = (180.0 + 20.) * Units.kg
    payload.mass_properties.center_of_gravity         = np.array([1.200,0.0,-0.112])
    net.payload                  = payload
    
    # Component 7 the Avionics
    avionics = SUAVE.Components.Energy.Peripherals.Avionics()
    avionics.power_draw = 0. * Units.watts
    avionics.mass_properties.mass = 10.0 * Units.kg
    avionics.mass_properties.center_of_gravity         = np.array([0.645,0.0,0.0])
    net.avionics        = avionics      

    # Component 8 the Battery
    #https://www.orbtronic.com/sony-vc7-18650-3500mah-li-ion-rechargeable-battery-green-flat-top
    bat = SUAVE.Components.Energy.Storages.Batteries.Constant_Mass.Lithium_Ion()

    bat.mass_properties.mass = 185.5 * Units.kg
    bat.mass_properties.center_of_gravity         = np.array([1.1,0.0,-0.1])
    casing_mass_fraction     = 0.1 # 10 percent casing mass fraction
    bat.system_voltage       = 120. * Units['volt']

    cell_voltage             = 3.6 * Units['volt']
    cell_Ah                  = 3.5 * Units['ampere'] / Units['hours']
    cell_mass                = 0.0471 * Units.kg
    cell_specific_energy     = cell_voltage * cell_Ah / cell_mass * Units.Wh/Units.kg
    #print('cell specific_energy: '+ str(cell_specific_energy) + ' Wh/kg')
    cell_resistance          = 23e-3
    
    bat.specific_energy      = cell_specific_energy * (1.-casing_mass_fraction) * Units.Wh/Units.kg
    #print('Bat. specific energy: ' + str(bat.specific_energy) + ' Wh/kg')

    number_cells         = np.round(bat.mass_properties.mass * (1.-casing_mass_fraction) / cell_mass)
    serial               = np.round(bat.system_voltage / cell_voltage)
    parallel             = np.round(number_cells / serial)

    #print('serial: ' + str(serial))
    #print('parallel: ' + str(parallel)) 

    bat.resistance           = 1/(parallel * 1/(serial * cell_resistance)) 
    #print('total battery resistance: ' + str(bat.resistance) + ' Ohm')
    #bat.resistance           = 0.003478 # resistance for 33s 115p Sony Murata Cells 
    
    initialize_from_mass(bat,bat.mass_properties.mass)
    net.battery              = bat
   
    #Component 9 the system logic controller and MPPT
    logic = SUAVE.Components.Energy.Distributors.Solar_Logic()
    logic.system_voltage  = 120.0
    logic.MPPT_efficiency = 0.965
    logic.mass_properties.mass = 3.0  * Units.kg
    
    net.solar_logic       = logic
    
    # add the solar network to the vehicle


    vehicle.append_component(net)
    
    # ------------------------------------------------------------------
    #   Vehicle Definition Complete
    # ------------------------------------------------------------------

    return vehicle

# ----------------------------------------------------------------------
#   Define the Configurations
# ---------------------------------------------------------------------

def configs_setup(vehicle):
    
    # ------------------------------------------------------------------
    #   Initialize Configurations
    # ------------------------------------------------------------------
    configs = SUAVE.Components.Configs.Config.Container()

    base_config = SUAVE.Components.Configs.Config(vehicle)
    base_config.tag = 'base'
    configs.append(base_config)

    return configs

def simple_sizing(configs):

    base = configs.base
    base.pull_base()

    # zero fuel weight
    base.mass_properties.max_zero_fuel = 0.9 * base.mass_properties.max_takeoff 

    # wing areas
    for wing in base.wings:
        wing.areas.wetted   = 2.0 * wing.areas.reference
        wing.areas.exposed  = 0.8 * wing.areas.wetted
        wing.areas.affected = 0.6 * wing.areas.wetted

    # diff the new data
    base.store_diff()

    return

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------

def mission_setup(analyses):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = SUAVE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'the_mission'

    #airport
    airport = SUAVE.Attributes.Airports.Airport()
    airport.altitude   =  0.0  * Units.ft
    airport.delta_isa  =  0.0
    airport.atmosphere = SUAVE.Attributes.Atmospheres.Earth.US_Standard_1976()

    mission.airport = airport    

    # unpack Segments module
    Segments = SUAVE.Analyses.Mission.Segments

    # base segment
    base_segment = Segments.Segment()

    # ------------------------------------------------------------------
    #   First Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1"
    segment.analyses.extend( analyses.base)  
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 7. * Units.deg    
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 3.0   * Units.km
    segment.air_speed      = 125.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to misison
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Climb Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg  
    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 190.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Climb Segment: constant Speed, Constant Rate
    # ------------------------------------------------------------------    

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg  
    segment.altitude_end = 10.668 * Units.km
    segment.air_speed    = 226.0  * Units['m/s']
    segment.climb_rate   = 3.0    * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------    
    #   Cruise Segment: Constant Speed, Constant Altitude
    # ------------------------------------------------------------------    

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"
    segment.analyses.extend( analyses.base)
    segment.air_speed  = 230.412 * Units['m/s']
    segment.distance   = 2490. * Units.nautical_miles

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg      
    segment.altitude_end = 8.0   * Units.km
    segment.air_speed    = 220.0 * Units['m/s']
    segment.descent_rate = 4.5   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 6.0   * Units.km
    segment.air_speed    = 195.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 4.0   * Units.km
    segment.air_speed    = 170.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fourth Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_4"
    segment.analyses.extend( analyses.base)
    segment.altitude_end = 2.0   * Units.km
    segment.air_speed    = 150.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fifth Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_5"
    segment.analyses.extend( analyses.base)
    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg       
    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 145.0 * Units['m/s']
    segment.descent_rate = 3.0   * Units['m/s']

    # append to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Mission definition complete    
    # ------------------------------------------------------------------

    return mission

def missions_setup(base_mission):

    # the mission container
    missions = SUAVE.Analyses.Mission.Mission.Container()

    # ------------------------------------------------------------------
    #   Base Mission
    # ------------------------------------------------------------------

    missions.base = base_mission

    return missions  

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,line_style='bo-'):

    axis_font = {'fontname':'Arial', 'size':'14'}    

    # ------------------------------------------------------------------
    #   Aerodynamics
    # ------------------------------------------------------------------


    fig = plt.figure("Aerodynamic Forces",figsize=(8,6))
    for segment in results.segments.values():

        time   = segment.conditions.frames.inertial.time[:,0] / Units.min
        Thrust = segment.conditions.frames.body.thrust_force_vector[:,0] / Units.lbf
        eta    = segment.conditions.propulsion.throttle[:,0]

        axes = fig.add_subplot(2,1,1)
        axes.plot( time , Thrust , line_style )
        axes.set_ylabel('Thrust (lbf)',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(2,1,2)
        axes.plot( time , eta , line_style )
        axes.set_xlabel('Time (min)',axis_font)
        axes.set_ylabel('Throttle',axis_font)
        axes.grid(True)	

        plt.savefig("B737_engine.pdf")
        plt.savefig("B737_engine.png")

    # ------------------------------------------------------------------
    #   Aerodynamics 2
    # ------------------------------------------------------------------
    fig = plt.figure("Aerodynamic Coefficients",figsize=(8,10))
    for segment in results.segments.values():

        time   = segment.conditions.frames.inertial.time[:,0] / Units.min
        CLift  = segment.conditions.aerodynamics.lift_coefficient[:,0]
        CDrag  = segment.conditions.aerodynamics.drag_coefficient[:,0]
        aoa = segment.conditions.aerodynamics.angle_of_attack[:,0] / Units.deg
        l_d = CLift/CDrag

        axes = fig.add_subplot(3,1,1)
        axes.plot( time , CLift , line_style )
        axes.set_ylabel('Lift Coefficient',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(3,1,2)
        axes.plot( time , l_d , line_style )
        axes.set_ylabel('L/D',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(3,1,3)
        axes.plot( time , aoa , 'ro-' )
        axes.set_xlabel('Time (min)',axis_font)
        axes.set_ylabel('AOA (deg)',axis_font)
        axes.grid(True)

        plt.savefig("B737_aero.pdf")
        plt.savefig("B737_aero.png")

    # ------------------------------------------------------------------
    #   Aerodynamics 2
    # ------------------------------------------------------------------
    fig = plt.figure("Drag Components",figsize=(8,10))
    axes = plt.gca()
    for i, segment in enumerate(results.segments.values()):

        time   = segment.conditions.frames.inertial.time[:,0] / Units.min
        drag_breakdown = segment.conditions.aerodynamics.drag_breakdown
        cdp = drag_breakdown.parasite.total[:,0]
        cdi = drag_breakdown.induced.total[:,0]
        cdc = drag_breakdown.compressible.total[:,0]
        cdm = drag_breakdown.miscellaneous.total[:,0]
        cd  = drag_breakdown.total[:,0]

        if line_style == 'bo-':
            axes.plot( time , cdp , 'ko-', label='CD parasite' )
            axes.plot( time , cdi , 'bo-', label='CD induced' )
            axes.plot( time , cdc , 'go-', label='CD compressibility' )
            axes.plot( time , cdm , 'yo-', label='CD miscellaneous' )
            axes.plot( time , cd  , 'ro-', label='CD total'   )
            if i == 0:
                axes.legend(loc='upper center')            
        else:
            axes.plot( time , cdp , line_style )
            axes.plot( time , cdi , line_style )
            axes.plot( time , cdc , line_style )
            axes.plot( time , cdm , line_style )
            axes.plot( time , cd  , line_style )            

    axes.set_xlabel('Time (min)')
    axes.set_ylabel('CD')
    axes.grid(True)
    plt.savefig("B737_drag.pdf")
    plt.savefig("B737_drag.png")

    # ------------------------------------------------------------------
    #   Altitude, sfc, vehicle weight
    # ------------------------------------------------------------------

    fig = plt.figure("Altitude_sfc_weight",figsize=(8,10))
    for segment in results.segments.values():

        time     = segment.conditions.frames.inertial.time[:,0] / Units.min
        aoa      = segment.conditions.aerodynamics.angle_of_attack[:,0] / Units.deg
        mass     = segment.conditions.weights.total_mass[:,0] / Units.lb
        altitude = segment.conditions.freestream.altitude[:,0] / Units.ft
        mdot     = segment.conditions.weights.vehicle_mass_rate[:,0]
        thrust   =  segment.conditions.frames.body.thrust_force_vector[:,0]
        sfc      = (mdot / Units.lb) / (thrust /Units.lbf) * Units.hr

        axes = fig.add_subplot(3,1,1)
        axes.plot( time , altitude , line_style )
        axes.set_ylabel('Altitude (ft)',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(3,1,3)
        axes.plot( time , sfc , line_style )
        axes.set_xlabel('Time (min)',axis_font)
        axes.set_ylabel('sfc (lb/lbf-hr)',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(3,1,2)
        axes.plot( time , mass , 'ro-' )
        axes.set_ylabel('Weight (lb)',axis_font)
        axes.grid(True)

        plt.savefig("B737_mission.pdf")
        plt.savefig("B737_mission.png")
        
    # ------------------------------------------------------------------
    #   Velocities
    # ------------------------------------------------------------------
    fig = plt.figure("Velocities",figsize=(8,10))
    for segment in results.segments.values():

        time     = segment.conditions.frames.inertial.time[:,0] / Units.min
        Lift     = -segment.conditions.frames.wind.lift_force_vector[:,2]
        Drag     = -segment.conditions.frames.wind.drag_force_vector[:,0] / Units.lbf
        Thrust   = segment.conditions.frames.body.thrust_force_vector[:,0] / Units.lb
        velocity = segment.conditions.freestream.velocity[:,0]
        pressure = segment.conditions.freestream.pressure[:,0]
        density  = segment.conditions.freestream.density[:,0]
        EAS      = velocity * np.sqrt(density/1.225)
        mach     = segment.conditions.freestream.mach_number[:,0]

        axes = fig.add_subplot(3,1,1)
        axes.plot( time , velocity / Units.kts, line_style )
        axes.set_ylabel('velocity (kts)',axis_font)
        axes.grid(True)

        axes = fig.add_subplot(3,1,2)
        axes.plot( time , EAS / Units.kts, line_style )
        axes.set_xlabel('Time (min)',axis_font)
        axes.set_ylabel('Equivalent Airspeed',axis_font)
        axes.grid(True)    
        
        axes = fig.add_subplot(3,1,3)
        axes.plot( time , mach , line_style )
        axes.set_xlabel('Time (min)',axis_font)
        axes.set_ylabel('Mach',axis_font)
        axes.grid(True)           
        
    return

if __name__ == '__main__': 
    main()    
    plt.show()
