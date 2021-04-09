# Jet_Aircraft_LJ45.py
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
ca.wingloadinglist_pa = np.arange(100, 1000, 5)* Units['force_kilogram/meter**2']
## Brief
# takeoff
ca.runway_elevation = 0*Units['meter']
ca.ground_run = 0.9*2180*Units['meters']  # takeoff distance given 2180
#climb
ca.climb_altitude = 0*Units['meter'] 
ca.climb_speed = 143*Units['knot'] # jenkinson
ca.climb_rate = 2500*Units['feet/minute'] # eurocontrol
#cruise
ca.cruise_altitude = 37000*Units['feet'] # jenkinson
ca.cruise_speed = 448*Units['knot'] # jenkinson
ca.cruise_thrust_fraction = 1.0
#turn
ca.turn_angle = 5*Units['degrees'] # should be enough for changing the heading
#ca.turn_load_factor = -1 #1.1547
ca.turn_altitude = 37000*Units['feet'] 
ca.turn_speed = 435*Units['knot'] #a slightly lower than cruise
ca.turn_specificenergy = 0
#ceiling
ca.service_ceiling = 39000*Units['feet']
ca.service_ceiling_climb_speed = 394*Units['knot'] # from flight manual 
#landing
ca.landing_roll = 0.6*1440*Units['meters'] # landing distance given 1440 


## Design
ca.aspect_ratio = 9.39  # jenkinson
ca.thickness_to_chord = 0.1192  # wikipedia reference paper
ca.oswald_factor = 0.8  # in paper stochastic drag polar 0.9220
# engine
ca.bypass_ratio = 6
ca.throttle_ratio = 1.07
# sweep
ca.sweep_leading_edge = 25*Units['degree']
#ca.sweep_max_thickness = 25*Units['degrees']
ca.sweep_quarter_chord = 25*Units['degree']

# high-lift 
ca.high_lift_type_takeoff = 3 #just a number 
ca.high_lift_type_landing = 3  #just a number 
ca.high_lift_type_clean = 0 
# weight
ca.climb_weight_fraction = 1.0
ca.cruise_weight_fraction = 70/77 # # should be reasonable from flight manual
ca.turn_weight_fraction = 70/77 # should be reasonable from flight manual
ca.ceiling_weight_fraction = 62/77 # flight manual 

## Performance
# take-off
ca.rolling_resistance = 0.05
ca.cd_takeoff = 0.078    # from paper stochastic drag polar
ca.cl_takeoff = 1.4     # random number
ca.cl_max_takeoff = 2.56    # jenkinson
# cruise, turn, ceiling
ca.cd_min_clean = 0.020  # in paper stochastic drag polar it was 0.023, other give a value around 0.018
#ca.cl_max_clean = 1
# landing
ca.cl_max_approach = 3.00  # jenkinson
# propeller
ca.propeller_efficiency_takeoff =  1
ca.propeller_efficiency_climb = 1
ca.propeller_efficiency_cruise = 1
ca.propeller_efficiency_turn = 1
ca.propeller_efficiency_servceil = 1 



#plot_thrust = ca.plot_constraint_diagram('thrust', 'US')
#plot_power = ca.plot_constraint_diagram('power', 'US')
a = ca.calculate_power_to_weight()
plot_power = ca.plot_constraint_diagram('thrust', 'si-kg')