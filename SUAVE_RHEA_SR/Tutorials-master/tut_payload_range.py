# tut_payload_range.py
#
# Created:  Aug 2014, SUAVE Team
# Modified: Apr 2016, T. Orra
#           Aug 2017, E. Botero

# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data 
from SUAVE.Methods.Propulsion.turbofan_sizing import turbofan_sizing
from SUAVE.Methods.Performance  import payload_range

import numpy as np
import pylab as plt

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    # define the problem
    configs, analyses = full_setup()

    configs.finalize()
    analyses.finalize()

    # weight analysis
    weights = analyses.configs.base.weights
    breakdown = weights.evaluate()

    # mission analysis
    mission = analyses.missions
    results = mission.evaluate()

    # run payload diagram
    config = configs.base
    cruise_segment_tag = "cruise"
    reserves = 1750.
    payload_range_results = payload_range(config,mission,cruise_segment_tag,reserves)

    # plot the results
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

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = mission

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

    # adjust analyses for configs

    # takeoff_analysis
    analyses.takeoff.aerodynamics.drag_coefficient_increment = 0.1000

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
    #weights = SUAVE.Analyses.Weights.Weights_Tube_Wing()
    weights = SUAVE.Analyses.Weights.Weights_FLOPS()
    weights.vehicle = vehicle

    # Define fuselage structural weight reduction factors
    weights.vehicle.settings = Data()
    weights.vehicle.settings.weight_reduction_factors = Data()
    weights.vehicle.settings.weight_reduction_factors.main_wing = 0.20          # 20% structure weight reduction in main wing 
    weights.vehicle.settings.weight_reduction_factors.fuselage  = 0.20          # 20% structure weight reduction in cabin and aft centerbody
    weights.vehicle.settings.weight_reduction_factors.empennage = 0.20          # 20% structure weight reduction in  HTP, VTP        
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.Fidelity_Zero()
    aerodynamics.geometry = vehicle
    
    # Define boundary layer control modifications (new transiton locations)
    aerodynamics.settings.fuselage_xt = 0.7                                 # Fuselage transition location
    aerodynamics.settings.wing_xt     = 0.8                                 # Wing transition location    
    aerodynamics.settings.prop_xt     = 0.7                                 # Propulstion system transition location 
 
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    #stability = SUAVE.Analyses.Stability.Fidelity_Zero()
    #stability.geometry = vehicle
    #analyses.append(stability)

    # ------------------------------------------------------------------
    #  Propulsion Analysis
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
    vehicle.tag = 'SE2A'

    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------
    # structure design conditions
    vehicle.type            = 'transport'
    vehicle.cruise_pressure = 26200.8 * Units.Pa
    vehicle.SL_pressure     = 101325.0 * Units.Pa
    vehicle.cruise_mach     = 0.78
    vehicle.max_cruise_mach = 0.84
    vehicle.range           = 2490. * Units.nautical_miles

    vehicle.movable_surface_area = 0.                                              # Total area of movable surfaces (flaps, slats, ailerons, elevator, rudder, etc)

    vehicle.fuel_tanks           = 1
    
    # mass properties
    vehicle.mass_properties.max_takeoff               = 57713.0 * Units.kg
    vehicle.mass_properties.operating_empty           = 32016.0 * Units.kg
    vehicle.mass_properties.takeoff                   = 57713.0 * Units.kg
    vehicle.mass_properties.max_zero_fuel             = 51597.0 * Units.kg
    vehicle.mass_properties.cargo                     = 0.0     * Units.kg
    vehicle.mass_properties.max_payload               = 18580.  * Units.kg
    vehicle.mass_properties.max_fuel                  = 10116.  * Units.kg
    vehicle.mass_properties.wing_cargo                = 0. * Units.kilogram
    vehicle.mass_properties.fuselage_cargo            = 0. * Units.kilogram
    vehicle.mass_properties.auxilary_fuel             = 0. * Units.kilogram

    # envelope properties
    vehicle.envelope.ultimate_load = 1.5
    vehicle.envelope.limit_load    = 1.0

    # basic parameters
    vehicle.reference_area         = 135 * Units['meters**2'] 
    vehicle.concept                ='TAW' # added a/c category for identifying TAW and BWB
    vehicle.crew                   = 2
    
    # 1-class Cabin Configuration
    vehicle.passengers             = 196
    vehicle.passengers_first       = 0
    vehicle.passengers_business    = 0
    vehicle.passengers_economy     = 196
    vehicle.seat_pitch_first       = 39 * Units.inches
    vehicle.seat_pitch_business    = 34 * Units.inches
    vehicle.seat_pitch_economy     = 31 * Units.inches
    vehicle.seats_per_row_first    = 4
    vehicle.seats_per_row_business = 6
    vehicle.seats_per_row_economy  = 6
    
    vehicle.systems.control        = "fully powered" 
    vehicle.systems.accessories    = "medium range"    

    # ------------------------------------------------------------------        
    #  Landing Gear
    # ------------------------------------------------------------------        
    # used for noise calculations
    landing_gear = SUAVE.Components.Landing_Gear.Landing_Gear()
    landing_gear.tag = "main_landing_gear"
    
    landing_gear.main_tire_diameter = 1.12000 * Units.meter
    landing_gear.nose_tire_diameter = 0.6858 * Units.meter
    landing_gear.main_strut_length  = 1.8 * Units.meter
    landing_gear.nose_strut_length  = 1.3 * Units.meter
    landing_gear.main_units  = 2                                    # number of main landing gear units
    landing_gear.nose_units  = 1                                    # number of nose landing gear
    landing_gear.main_wheels = 2                                    # number of wheels on the main landing gear
    landing_gear.nose_wheels = 2                                    # number of wheels on the nose landing gear      
    vehicle.landing_gear     = landing_gear
    
    # ------------------------------------------------------------------
    #   Main Wing
    # ------------------------------------------------------------------

    wing = SUAVE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'
    
    wing.aspect_ratio            = 16
    wing.sweeps.quarter_chord    = -17 * Units.degrees
    wing.dihedral                = 3.0 * Units.degrees
    wing.thickness_to_chord      = 0.1
    wing.taper                   = 0.25
    wing.span_efficiency         = 0.9
    wing.spans.projected         = 46.3 * Units.meter
    wing.chords.root             = 4.65 * Units.meter
    wing.chords.tip              = 1.16 * Units.meter
    wing.chords.mean_aerodynamic = 3.25 * Units.meter
    wing.areas.reference         = 135  * Units['meters**2']
    wing.areas.flap              = 0.24 * wing.areas.reference
    wing.twists.root             = 2.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees
    wing.origin                  = [20.5,0,-1.27]              # meters
    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = True
    wing.folding                 = True
    wing.folding_penalty         = 0.21                         # Weight penalty of the folding mechanism (% of the wing weight)
    wing.FSTRT                   = 0.0                          # Wing strut bracing factor (0 - no bracing to 1 - full bracing)
    wing.FAERT                   = 1.0                          # Aeroelastic tailoring factor (0 - no tailoring tp 1 - max tailoring)
    wing.PCTL                    = 1.0                          # Load percent distribution 
    wing.dynamic_pressure_ratio  = 1.0 

    # add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #  Horizontal Stabilizer
    # ------------------------------------------------------------------

    wing = SUAVE.Components.Wings.Wing()
    wing.tag = 'horizontal_stabilizer'

    wing.aspect_ratio            = 6.53    
    wing.sweeps.quarter_chord    = 15 * Units.degrees
    wing.thickness_to_chord      = 0.09
    wing.taper                   = 0.4
    wing.span_efficiency         = 0.9
    wing.spans.projected         = 14.9 * Units.meter
    wing.chords.root             = 3.5  * Units.meter
    wing.chords.tip              = 1.3 * Units.meter
    wing.chords.mean_aerodynamic = 2.43  * Units.meter
    wing.areas.reference         = 34.14   * Units['meters**2']   
    wing.areas.wetted            = 2.0  * wing.areas.reference
    wing.areas.exposed           = 0.8  * wing.areas.wetted
    wing.areas.affected          = 0.6  * wing.areas.wetted
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [32.83,0,1.14] # meters
    wing.vertical                = False 
    wing.symmetric               = True
    wing.dynamic_pressure_ratio  = 0.9
    
    # add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------

    wing = SUAVE.Components.Wings.Wing()
    wing.tag = 'vertical_stabilizer'

    wing.aspect_ratio            = 1.57
    wing.sweeps.quarter_chord    = 15 * Units.degrees
    wing.thickness_to_chord      = 0.09
    wing.taper                   = 0.6
    wing.span_efficiency         = 0.9
    wing.spans.projected         = 7.5 * Units.meter
    wing.chords.root             = 5.85  * Units.meter
    wing.chords.tip              = 3.5  * Units.meter
    wing.chords.mean_aerodynamic = 4.78   * Units.meter
    wing.areas.reference         = 35.8 * Units['meters**2']  
    wing.areas.wetted            = 2.0  * wing.areas.reference
    wing.areas.exposed           = 0.8  * wing.areas.wetted
    wing.areas.affected          = 0.6  * wing.areas.wetted
    wing.twists.root             = 0.0 * Units.degrees
    wing.twists.tip              = 0.0 * Units.degrees  
    wing.origin                  = [32,0,1.54] # meters
    wing.vertical                = True 
    wing.symmetric               = False
    wing.t_tail                  = False
    wing.dynamic_pressure_ratio  = 1.0

    # add to vehicle
    vehicle.append_component(wing)

    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------

    fuselage = SUAVE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage'
    
    fuselage.number_coach_seats    = vehicle.passengers
    fuselage.seats_abreast         = 6
    fuselage.number_of_aisles      = 1
    fuselage.seat_pitch            = 0.86   * Units.meter
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
    fuselage.effective_diameter    = 3.74  * Units.meter
    fuselage.areas.top_projected   = 142.1948 * Units['meters**2'] 
    fuselage.areas.side_projected  = 142.1948 * Units['meters**2']
    fuselage.areas.cabin_floor     = 100.00   * Units['meters**2']
    fuselage.areas.wetted          = 446.718  * Units['meters**2'] 
    fuselage.areas.front_projected = 12.57    * Units['meters**2'] 
    fuselage.differential_pressure = 5.0e4    * Units.pascal             # Maximum differential pressure
    
    fuselage.heights.at_quarter_length          = 3.74 * Units.meter
    fuselage.heights.at_three_quarters_length   = 3.65 * Units.meter
    fuselage.heights.at_wing_root_quarter_chord = 3.74 * Units.meter    

    # add to vehicle
    vehicle.append_component(fuselage)

    #-------------------------------------------------------------------
    #   Paint
    #-------------------------------------------------------------------
    vehicle.paint = Data()
    vehicle.paint.area_density = 0.5     # Units.kg/Units['meters**2']

    
    # ------------------------------------------------------------------
    #  Turbofan Network
    # ------------------------------------------------------------------


    # initialize the gas turbine network
    gt_engine                   = SUAVE.Components.Energy.Networks.Turbofan()
    gt_engine.tag               = 'turbofan'

    # setup
    gt_engine.number_of_engines        = 3
    gt_engine.wing_mounted_engines     = 2
    gt_engine.fuselage_mounted_engines = 1
    gt_engine.engine_position          = 'on wings'                         # position of the engines ('on wings' if any are on wings, 'otherwise' - in none are on the wings)
    gt_engine.bypass_ratio             = 5.4
    gt_engine.engine_length            = 2.71 * Units.meter
    gt_engine.nacelle_diameter         = 2.05 * Units.meter
    gt_engine.origin                   = [[13.72, 4.86,-1.9],[13.72, -4.86,-1.9]] # meters

    # engine weights
    gt_engine.base_engine_weight   = 696 * Units.kilogram              # weight of a selected engine
    gt_engine.nozzle_included      = 0                               # indicator if the nozzle is included in weights (0 - no, 1 - yes) 
    gt_engine.nozzle_weight        = 0 * Units.kilogram
    gt_engine.inlet_included       = 0                               # indicator if the inlet is included in weights (0 - no, 1 - yes) 
    gt_engine.inlet_weight         = 0 * Units.kilogram
    gt_engine.miscellaneous_weight = 0 * Units.kilogram              # miscellaneous weight components weights

    # compute engine areas
    gt_engine.areas.wetted      = 1.1*np.pi*gt_engine.nacelle_diameter*gt_engine.engine_length

    # set the working fluid for the network
    gt_engine.working_fluid = SUAVE.Attributes.Gases.Air()

    # Component 1 : ram,  to convert freestream static to stagnation quantities
    ram = SUAVE.Components.Energy.Converters.Ram()
    ram.tag = 'ram'

    # add ram to the network
    gt_engine.ram = ram

    # Component 2 : inlet nozzle
    inlet_nozzle = SUAVE.Components.Energy.Converters.Compression_Nozzle()
    inlet_nozzle.tag = 'inlet nozzle'

    inlet_nozzle.polytropic_efficiency = 0.98
    inlet_nozzle.pressure_ratio        = 0.99

    # add inlet nozzle to the network
    gt_engine.inlet_nozzle = inlet_nozzle

    # Component 3 :low pressure compressor
    low_pressure_compressor = SUAVE.Components.Energy.Converters.Compressor()
    low_pressure_compressor.tag = 'lpc'

    low_pressure_compressor.polytropic_efficiency = 0.91
    low_pressure_compressor.pressure_ratio        = 1.9

    # add low pressure compressor to the network
    gt_engine.low_pressure_compressor = low_pressure_compressor

    # Component 4 :high pressure compressor
    high_pressure_compressor = SUAVE.Components.Energy.Converters.Compressor()
    high_pressure_compressor.tag = 'hpc'

    high_pressure_compressor.polytropic_efficiency = 0.91
    high_pressure_compressor.pressure_ratio        = 10.0

    # add the high pressure compressor to the network
    gt_engine.high_pressure_compressor = high_pressure_compressor

    # Component 5 :low pressure turbine
    low_pressure_turbine = SUAVE.Components.Energy.Converters.Turbine()
    low_pressure_turbine.tag='lpt'

    low_pressure_turbine.mechanical_efficiency = 0.99
    low_pressure_turbine.polytropic_efficiency = 0.99

    # add low pressure turbine to the network
    gt_engine.low_pressure_turbine = low_pressure_turbine

    # Component 5 :high pressure turbine
    high_pressure_turbine = SUAVE.Components.Energy.Converters.Turbine()
    high_pressure_turbine.tag='hpt'

    high_pressure_turbine.mechanical_efficiency = 0.99
    high_pressure_turbine.polytropic_efficiency = 0.99

    # add the high pressure turbine to the network
    gt_engine.high_pressure_turbine = high_pressure_turbine

    # Component 6 :combustor
    combustor = SUAVE.Components.Energy.Converters.Combustor()
    combustor.tag = 'comb'
    combustor.efficiency                = 0.99
    combustor.alphac                    = 1.0
    combustor.turbine_inlet_temperature = 1500
    combustor.pressure_ratio            = 0.95
    combustor.fuel_data                 = SUAVE.Attributes.Propellants.Jet_A()

    # add the combustor to the network
    gt_engine.combustor = combustor

    # Component 7 :core nozzle
    core_nozzle = SUAVE.Components.Energy.Converters.Expansion_Nozzle()
    core_nozzle.tag = 'core nozzle'

    core_nozzle.polytropic_efficiency = 0.95
    core_nozzle.pressure_ratio        = 0.99

    # add the core nozzle to the network
    gt_engine.core_nozzle = core_nozzle

    # Component 8 :fan nozzle
    fan_nozzle = SUAVE.Components.Energy.Converters.Expansion_Nozzle()
    fan_nozzle.tag = 'fan nozzle'

    fan_nozzle.polytropic_efficiency = 0.95
    fan_nozzle.pressure_ratio        = 0.98

    # add the fan nozzle to the network
    gt_engine.fan_nozzle = fan_nozzle

    # Component 9 : fan
    fan = SUAVE.Components.Energy.Converters.Fan()
    fan.tag = 'fan'

    fan.polytropic_efficiency = 0.93
    fan.pressure_ratio        = 1.7

    # add the fan to the network
    gt_engine.fan = fan

    # Component 10 : thrust reverser (to compute landing performance and weights)
    gt_engine.thrust_reverser = True

    # Component 11 : thrust (to compute the thrust)
    thrust = SUAVE.Components.Energy.Processes.Thrust()
    thrust.tag ='compute_thrust'

    # total design thrust (includes all the engines)
    thrust.total_design             = 3*25000. * Units.N 
    thrust.baseline                 = 3*33600. * Units.N           # Baseline engine thrust at SL (if no engine was selected, use 0)

    # design sizing conditions
    altitude      = 33000.0*Units.ft
    mach_number   = 0.78
    isa_deviation = 0.

    # add thrust to the network
    gt_engine.thrust = thrust

    # size the turbofan
    turbofan_sizing(gt_engine,mach_number,altitude)

    # add  gas turbine network gt_engine to the vehicle
    vehicle.append_component(gt_engine)

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

    # ------------------------------------------------------------------
    #   Cruise Configuration
    # ------------------------------------------------------------------

    config = SUAVE.Components.Configs.Config(base_config)
    config.tag = 'cruise'

    configs.append(config)

    # ------------------------------------------------------------------
    #   Takeoff Configuration
    # ------------------------------------------------------------------

    config = SUAVE.Components.Configs.Config(base_config)
    config.tag = 'takeoff'

    config.wings['main_wing'].flaps.angle = 20. * Units.deg
    config.wings['main_wing'].slats.angle = 25. * Units.deg

    config.V2_VS_ratio = 1.21
    config.maximum_lift_coefficient = 2.

    configs.append(config)

    # ------------------------------------------------------------------
    #   Landing Configuration
    # ------------------------------------------------------------------

    config = SUAVE.Components.Configs.Config(base_config)
    config.tag = 'landing'

    config.wings['main_wing'].flaps_angle = 30. * Units.deg
    config.wings['main_wing'].slats_angle = 25. * Units.deg

    config.Vref_VS_ratio = 1.23
    config.maximum_lift_coefficient = 2.

    configs.append(config)

    return configs

# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------
def mission_setup(analyses):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    mission = SUAVE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'embraer_e190ar test mission'

    # atmospheric model
    atmosphere = SUAVE.Attributes.Atmospheres.Earth.US_Standard_1976()
    planet = SUAVE.Attributes.Planets.Earth()

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
    #   First Climb Segment
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_1"

    # connect vehicle configuration
    segment.analyses.extend( analyses.takeoff )

    # define segment attributes
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 3.0   * Units.km
    segment.air_speed      = 125.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to misison
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Climb Segment
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_2"

    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )

    # segment attributes
    segment.altitude_end   = 8.0   * Units.km
    segment.air_speed      = 190.0 * Units['m/s']
    segment.climb_rate     = 6.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Third Climb Segment
    # ------------------------------------------------------------------

    segment = Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb_3"

    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )

    # segment attributes
    segment.altitude_end = 10.668 * Units.km
    segment.air_speed    = 226.0  * Units['m/s']
    segment.climb_rate   = 3.0    * Units['m/s']

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Cruise Segment
    # ------------------------------------------------------------------

    segment = Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"

    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )

    # segment attributes
    segment.atmosphere = atmosphere
    segment.planet     = planet

    segment.air_speed  = 230.412 * Units['m/s']
    segment.distance   = 1750. * Units.nautical_miles

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   First Descent Segment
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_1"

    # connect vehicle configuration
    segment.analyses.extend( analyses.cruise )

    # segment attributes
    segment.altitude_end = 8.0   * Units.km
    segment.air_speed    = 220.0 * Units['m/s']
    segment.descent_rate = 4.5   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Second Descent Segment
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_2"

    # connect vehicle configuration
    segment.analyses.extend( analyses.landing )

    # segment attributes
    segment.altitude_end = 6.0   * Units.km
    segment.air_speed    = 195.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)


    # ------------------------------------------------------------------
    #   Third Descent Segment
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_3"

    # connect vehicle configuration
    segment.analyses.extend( analyses.landing )

    # segment attributes
    segment.altitude_end = 4.0   * Units.km
    segment.air_speed    = 170.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']

    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fourth Descent Segment
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_4"

    segment.analyses.extend( analyses.landing )

    segment.altitude_end = 2.0   * Units.km
    segment.air_speed    = 150.0 * Units['m/s']
    segment.descent_rate = 5.0   * Units['m/s']


    # add to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Fifth Descent Segment
    # ------------------------------------------------------------------

    segment = Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent_5"

    segment.analyses.extend( analyses.landing )

    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 145.0 * Units['m/s']
    segment.descent_rate = 3.0   * Units['m/s']

    # append to mission
    mission.append_segment(segment)

    # ------------------------------------------------------------------
    #   Mission definition complete
    # ------------------------------------------------------------------

    return mission

# ----------------------------------------------------------------------
#   Plot Mission
# ----------------------------------------------------------------------

def plot_mission(results,line_style='bo-'):

    # ------------------------------------------------------------------
    #   Throttle
    # ------------------------------------------------------------------
    plt.figure("Throttle History")
    axes = plt.gca()
    for i in range(len(results.segments)):
        time = results.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        eta  = results.segments[i].conditions.propulsion.throttle[:,0]
        axes.plot(time, eta, line_style)
    axes.set_xlabel('Time (mins)')
    axes.set_ylabel('Throttle')
    axes.grid(True)

    # ------------------------------------------------------------------
    #   Angle of Attack
    # ------------------------------------------------------------------

    plt.figure("Angle of Attack History")
    axes = plt.gca()
    for i in range(len(results.segments)):
        time = results.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        aoa = results.segments[i].conditions.aerodynamics.angle_of_attack[:,0] / Units.deg
        axes.plot(time, aoa, line_style)
    axes.set_xlabel('Time (mins)')
    axes.set_ylabel('Angle of Attack (deg)')
    axes.grid(True)

    # ------------------------------------------------------------------
    #   Fuel Burn Rate
    # ------------------------------------------------------------------
    plt.figure("Fuel Burn Rate")
    axes = plt.gca()
    for i in range(len(results.segments)):
        time = results.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        mdot = results.segments[i].conditions.weights.vehicle_mass_rate[:,0]
        axes.plot(time, mdot, line_style)
    axes.set_xlabel('Time (mins)')
    axes.set_ylabel('Fuel Burn Rate (kg/s)')
    axes.grid(True)

    # ------------------------------------------------------------------
    #   Altitude
    # ------------------------------------------------------------------
    plt.figure("Altitude")
    axes = plt.gca()
    for i in range(len(results.segments)):
        time     = results.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        altitude = results.segments[i].conditions.freestream.altitude[:,0] / Units.km
        axes.plot(time, altitude, line_style)
    axes.set_xlabel('Time (mins)')
    axes.set_ylabel('Altitude (km)')
    axes.grid(True)

    # ------------------------------------------------------------------
    #   Vehicle Mass
    # ------------------------------------------------------------------
    plt.figure("Vehicle Mass")
    axes = plt.gca()
    for i in range(len(results.segments)):
        time = results.segments[i].conditions.frames.inertial.time[:,0] / Units.min
        mass = results.segments[i].conditions.weights.total_mass[:,0]
        axes.plot(time, mass, line_style)
    axes.set_xlabel('Time (mins)')
    axes.set_ylabel('Vehicle Mass (kg)')
    axes.grid(True)

    # ------------------------------------------------------------------
    #   Aerodynamics
    # ------------------------------------------------------------------
    fig = plt.figure("Aerodynamic Forces")
    for segment in results.segments.values():

        time   = segment.conditions.frames.inertial.time[:,0] / Units.min
        Lift   = -segment.conditions.frames.wind.lift_force_vector[:,2]
        Drag   = -segment.conditions.frames.wind.drag_force_vector[:,0]
        Thrust = segment.conditions.frames.body.thrust_force_vector[:,0]

        axes = fig.add_subplot(3,1,1)
        axes.plot( time , Lift , line_style )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('Lift (N)')
        axes.grid(True)

        axes = fig.add_subplot(3,1,2)
        axes.plot( time , Drag , line_style )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('Drag (N)')
        axes.grid(True)

        axes = fig.add_subplot(3,1,3)
        axes.plot( time , Thrust , line_style )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('Thrust (N)')
        axes.grid(True)

    # ------------------------------------------------------------------
    #   Aerodynamics 1
    # ------------------------------------------------------------------
    fig = plt.figure("Aerodynamic Coefficients")
    for segment in results.segments.values():

        time   = segment.conditions.frames.inertial.time[:,0] / Units.min
        CLift  = segment.conditions.aerodynamics.lift_coefficient[:,0]
        CDrag  = segment.conditions.aerodynamics.drag_coefficient[:,0]
        Drag   = -segment.conditions.frames.wind.drag_force_vector[:,0]
        Thrust = segment.conditions.frames.body.thrust_force_vector[:,0]

        axes = fig.add_subplot(3,1,1)
        axes.plot( time , CLift , line_style )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('CL')
        axes.grid(True)

        axes = fig.add_subplot(3,1,2)
        axes.plot( time , CDrag , line_style )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('CD')
        axes.grid(True)

        axes = fig.add_subplot(3,1,3)
        axes.plot( time , Drag   , line_style )
        axes.plot( time , Thrust , 'ro-' )
        axes.set_xlabel('Time (min)')
        axes.set_ylabel('Drag and Thrust (N)')
        axes.grid(True)


    # ------------------------------------------------------------------
    #   Aerodynamics 2
    # ------------------------------------------------------------------
    fig = plt.figure("Drag Components")
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
            axes.plot( time , cdp , 'ko-', label='CD_P' )
            axes.plot( time , cdi , 'bo-', label='CD_I' )
            axes.plot( time , cdc , 'go-', label='CD_C' )
            axes.plot( time , cdm , 'yo-', label='CD_M' )
            axes.plot( time , cd  , 'ro-', label='CD'   )
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

    return


if __name__ == '__main__':
    main()
    plt.show()