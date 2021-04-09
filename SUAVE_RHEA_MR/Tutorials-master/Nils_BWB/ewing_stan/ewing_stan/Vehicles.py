# Procedure.py
# 
# Created:  Apr 2019, N. Kleemann
# Modified: July 2019, N. Kleemann


#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data
import numpy as np
from SUAVE.Components.Energy.Networks.Solar import Solar
from SUAVE.Methods.Propulsion import propeller_design
from SUAVE.Methods.Power.Battery.Sizing import initialize_from_energy_and_power, initialize_from_mass



# ----------------------------------------------------------------------
#   Build the Vehicle
# ----------------------------------------------------------------------

def vehicle_setup():
    
    # ------------------------------------------------------------------
    #   Initialize the Vehicle
    # ------------------------------------------------------------------    
    
    vehicle = SUAVE.Vehicle()
    vehicle.tag = 'E-WingSolar'
    
    # ------------------------------------------------------------------
    #   Vehicle-level Properties
    # ------------------------------------------------------------------    
    # mass properties
    vehicle.mass_properties.takeoff             = 580. * Units.kg
    #vehicle.mass_properties.operating_empty    = 385. * Units.kg
    vehicle.mass_properties.max_takeoff         = 600. * Units.kg
    vehicle.mass_properties.payload             = 200. * Units.kg 

    # basic parameters
    vehicle.reference_area                      = 16.5 * Units["m^2"]
    #vehicle.envelope.ultimate_load             = 10.0
    #vehicle.envelope.limit_load                = 6.0
    #vehicle.envelope.maximum_dynamic_pressure  = 0.5*1.225*(100.**2.) #Max q

    vehicle.number_passengers                   = 2
    
    # ------------------------------------------------------------------        
    #   materials definition
    # ------------------------------------------------------------------ 

    bid_carbon      = SUAVE.Attributes.Solids.Bidirectional_Carbon_Fiber()
    bid_glass       = SUAVE.Attributes.Solids.Bidirectional_Glass_Fiber()

    ud_carbon       = SUAVE.Attributes.Solids.Unidirectional_Carbon_Fiber()
    ud_glass        = SUAVE.Attributes.Solids.Unidirectional_Glass_Fiber()

    core            = SUAVE.Attributes.Solids.Carbon_Fiber_Honeycomb()

    epoxy           = SUAVE.Attributes.Solids.Epoxy()
    paint           = SUAVE.Attributes.Solids.Paint()
    acrylic         = SUAVE.Attributes.Solids.Acrylic()

    aluminium_rib   = SUAVE.Attributes.Solids.Aluminum_Rib()
    glass_rib       = SUAVE.Attributes.Solids.Glass_Fiber_Rib()
    steel           = SUAVE.Attributes.Solids.Steel()

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


    airfoil = SUAVE.Components.Wings.Airfoils.Airfoil()
    airfoil.coordinate_file = 'NACA65_410.dat'

    wing = SUAVE.Components.Wings.Main_Wing()
    wing.tag = 'main_wing'

    wing.areas.reference         = vehicle.reference_area

    wing.spans.projected 		 = 11. * Units.meter
    wing.aspect_ratio            = (wing.spans.projected**2.) / wing.areas.reference

    wing.chords.root             = 3.4 * Units.meter
    wing.chords.seat             = 3.2 * Units.meter
    wing.chords.kink             = 1.55 * Units.meter
    wing.chords.tip              = 0.4 * Units.meter
    wing.chords.mid_wing         = (wing.chords.kink + wing.chords.tip) / 2. 

    spanwise_locations_root        = 0.0
    spanwise_locations_seat        = 0.6 * Units.meter / wing.spans.projected * 2
    spanwise_locations_kink        = 1.5 * Units.meter / wing.spans.projected * 2
    spanwise_locations_tip         = 1.0
    spanwise_locations_mid_wing    = (spanwise_locations_kink + spanwise_locations_tip) / 2

    wing.areas.reference         = ((wing.chords.root + wing.chords.seat) * spanwise_locations_seat +
                                    (wing.chords.seat + wing.chords.kink) * (spanwise_locations_kink - spanwise_locations_seat) +
                                    (wing.chords.tip + wing.chords.kink) * (spanwise_locations_tip - spanwise_locations_kink)) * wing.spans.projected / 2.

    print("wing areas reference: " + str(wing.areas.reference))

    vehicle.reference_area       = wing.areas.reference

    wing.taper                   = wing.chords.tip / wing.chords.root
    wing.chords.mean_aerodynamic = wing.areas.reference / wing.spans.projected
    #wing.sweeps.quarter_chord    = 25. * Units.degrees

    wing.twists.root             = -2.5 * Units.degrees
    wing.twists.seat             = -2.0 * Units.degrees
    wing.twists.kink             = 0.0 * Units.degrees
    wing.twists.tip              = -7.0 * Units.degrees
    wing.twists.mid_wing         = (wing.twists.tip - wing.twists.kink) / 4.
    
    wing.thickness_to_chord      = 0.18

    wing.origin                  = [0.,0.,0]
    wing.aerodynamic_center      = [0,0,0] 

    wing.vertical                = False
    wing.symmetric               = True
    wing.high_lift               = False

    wing.dynamic_pressure_ratio  = 1.0

    wing.transition_x_upper      = 0.5
    wing.transition_x_lower      = 0.8


    wing.number_ribs             = 6.
    wing.winglet_fraction        = 0.0
    wing.motor_spanwise_locations = [0.]


    wing.materials.skin_materials.torsion_carrier   = bid_glass
    wing.materials.skin_materials.core              = core
    wing.materials.skin_materials.adhesive          = epoxy
    wing.materials.skin_materials.covering          = paint

    wing.materials.flap_materials.bending_carrier   = ud_carbon

    wing.materials.rib_materials.structural         = glass_rib

    wing.materials.spar_materials.shear_carrier     = bid_glass
    wing.materials.spar_materials.bending_carrier   = ud_carbon



    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'root'
    segment.percent_span_location = spanwise_locations_root
    segment.twist                 = wing.twists.root
    segment.root_chord_percent    = 1.0
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 8.0 * Units.degrees
    segment.thickness_to_chord    = 0.32
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)    
    
    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'seat'
    segment.percent_span_location = spanwise_locations_seat
    segment.twist                 = wing.twists.seat
    segment.root_chord_percent    = wing.chords.seat / wing.chords.root
    segment.dihedral_outboard     = 1.   * Units.degrees
    segment.sweeps.quarter_chord  = 7.0 * Units.degrees
    segment.thickness_to_chord    = 0.335
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment)  


    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'kink'
    segment.percent_span_location = spanwise_locations_kink
    segment.twist                 = wing.twists.kink
    segment.root_chord_percent    = wing.chords.kink / wing.chords.root
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 25. * Units.degrees    
    segment.thickness_to_chord    = 0.18
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment) 

    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'mid_wing'
    segment.percent_span_location = spanwise_locations_mid_wing
    segment.twist                 = wing.twists.mid_wing
    segment.root_chord_percent    = wing.chords.mid_wing / wing.chords.root
    segment.dihedral_outboard     = 1. * Units.degrees
    segment.sweeps.quarter_chord  = 25. * Units.degrees
    segment.thickness_to_chord    = 0.16
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment) 

    
    segment = SUAVE.Components.Wings.Segment()
    segment.tag                   = 'tip'
    segment.percent_span_location = spanwise_locations_tip
    segment.twist                 = wing.twists.tip
    segment.root_chord_percent    = wing.chords.tip / wing.chords.root
    segment.sweeps.quarter_chord  = 25. * Units.degrees
    segment.dihedral_outboard = 1. * Units.degrees
    segment.append_airfoil(airfoil)
    wing.Segments.append(segment) 



    # ------------------------------------------------------------------
    #   Control Surfaces "Flaps"
    # ------------------------------------------------------------------
    wing.flaps.chord      =  0.25   # 30% of the chord
    wing.flaps.span_start =  0.7   # 10% of the span
    wing.flaps.span_end   =  1.0
    #wing.flaps.type       = 'double_slotted'

    # add to vehicle
    vehicle.append_component(wing)



    # ------------------------------------------------------------------
    #   Vertical Stabilizer
    # ------------------------------------------------------------------
    '''
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


    wing.materials.skin_materials.torsion_carrier   = bid_glass
    wing.materials.skin_materials.core              = core
    wing.materials.skin_materials.adhesive          = epoxy
    wing.materials.skin_materials.covering          = paint

    wing.materials.flap_materials.bending_carrier   = ud_carbon

    wing.materials.rib_materials.structural         = glass_rib

    wing.materials.spar_materials.shear_carrier     = bid_glass
    wing.materials.spar_materials.bending_carrier   = ud_carbon


    wing.transition_x_upper      = 0.5
    wing.transition_x_lower      = 0.5
        
    # add to vehicle
    vehicle.append_component(wing)
    '''


    # ------------------------------------------------------------------
    #  Fuselage
    # ------------------------------------------------------------------
    
    fuselage = SUAVE.Components.Fuselages.Fuselage()
    fuselage.tag = 'fuselage_bwb'

    fuselage.lengths.total         = 3.4 * Units.meter
    fuselage.width                 = 2.0  * Units.meter
    fuselage.heights.maximum       = 0.335 * 3.2 * Units.meter
    fuselage.areas.wetted          = 0.0  * Units['meters**2'] 
    fuselage.areas.front_projected = 0.0    * Units['meters**2']


    fuselage.materials.keel_materials.root_bending_moment_carrier = ud_glass
    fuselage.materials.keel_materials.shear_carrier = bid_glass
    fuselage.materials.keel_materials.bearing_carrier = ud_glass
    fuselage.materials.bolt_materials.landing_pad_bolt = steel
    fuselage.materials.canopy_materials.canopy = acrylic
    fuselage.materials.skin_materials.skin = bid_glass


    # add to vehicle
    vehicle.append_component(fuselage)

   
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


    prop.materials = SUAVE.Core.Container()
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
    prop.materials.skin_materials.cover                 = epoxy


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

    # ----------------------------------------------------------------------
    #   Center of Gravity and Stability
    # ---------------------------------------------------------------------
    #vehicle.mass_properties.mass                      = 600.0
    #vehicle.mass_properties.volume                    = 5.0
    #vehicle.mass_properties.center_of_gravity         = np.array([1.385,0.0,-0.148])
        
    #vehicle.mass_properties.moments_of_inertia.center = vehicle.center_of_gravity
    #vehicle.mass_properties.moments_of_inertia.tensor = np.array([[1239,0.0,7],[0.0,259,0.2],[7,0.2,1427]])

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
    
    return configs
