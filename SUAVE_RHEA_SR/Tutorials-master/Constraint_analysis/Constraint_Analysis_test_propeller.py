# Constraint_Analysis_test.py
# 
# Created:  Aug 2020, Radimir Yanev


# ----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import SUAVE.Analyses.Constraint_Analysis as suaveCA

from SUAVE.Core import Units


ca = suaveCA.Constraint_Analysis()
ca.wingloadinglist_pa = np.arange(5, 125, 5)* Units['force_pound/feet**2']
## Brief
# takeoff
ca.runway_elevation = 0*Units['meter']
ca.ground_run = 3800*Units['feet'] #1158.24 #3800 ft
#climb
ca.climb_altitude = 0*Units['meter']
ca.climb_speed = 255*Units['knot']
ca.climb_rate = 1800 *Units['feet/minute']
#cruise
ca.cruise_altitude = 20000*Units['feet']
ca.cruise_speed = 255*Units['knot']
ca.cruise_thrust_fraction = 1.0
#turn
ca.turn_angle = 30*Units['degrees']
#ca.turn_load_factor = -1 #1.1547
ca.turn_altitude = 20000*Units['feet'] #6096
ca.turn_speed = 255*Units['knot'] # 255*Units['knots']
ca.turn_specificenergy = 0
#ceiling
ca.service_ceiling = 24000*Units['feet']
ca.service_ceiling_climb_speed = 255*Units['knot']
#landing
ca.landing_roll = 2500*Units['feet'] # 2100*0.3048


## Design
ca.aspect_ratio = 14
ca.thickness_to_chord = 0.15
ca.oswald_factor = 0.844830315
# engine
ca.bypass_ratio = -1#-1
ca.throttle_ratio = 1.07
# sweep
ca.sweep_leading_edge = 0*Units['degree']
#ca.sweep_max_thickness = 0*Units['degrees']
ca.sweep_quarter_chord = 0

# high-lift
ca.high_lift_type_takeoff = 3
ca.high_lift_type_landing = 6
ca.high_lift_type_clean = 0
# weight
ca.climb_weight_fraction = 1.0
ca.cruise_weight_fraction = 1.0
ca.turn_weight_fraction = 1.0
ca.ceiling_weight_fraction = 1.0

## Performance
# take-off
ca.rolling_resistance = 0.05
ca.cd_takeoff = 0.12
ca.cl_takeoff = 1.16
ca.cl_max_takeoff = 2.2       
# cruise, turn, ceiling
ca.cd_min_clean = 0.0134
ca.cl_max_clean = 1.8999
# landing
ca.cl_max_approach = 1.8999       
# propeller
ca.propeller_efficiency_takeoff =  0.54
ca.propeller_efficiency_climb = 0.85
ca.propeller_efficiency_cruise = 0.85
ca.propeller_efficiency_turn = 0.85
ca.propeller_efficiency_servceil = 0.85



#plot_thrust = ca.plot_constraint_diagram('thrust', 'US')
#plot_power = ca.plot_constraint_diagram('power', 'US')
a = ca.calculate_power_to_weight()
plot_power = ca.plot_constraint_diagram('power', 'us')
