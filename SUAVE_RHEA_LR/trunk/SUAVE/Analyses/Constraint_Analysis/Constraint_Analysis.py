## @ingroup Analyses-Constraint_Analysis
# Constraint_Analysis.py
#
# Created: Aug 2020, R.Y. Yanev 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
import numpy as np
import matplotlib.pyplot as plt
from SUAVE.Core import Units

from SUAVE.Analyses.Constraint_Analysis.ADRpy import constraintanalysis as ca
from SUAVE.Analyses.Constraint_Analysis.ADRpy import atmospheres as at
# ----------------------------------------------------------------------
#  Analysis
# ----------------------------------------------------------------------
## @ingroup Analyses-Constraint_Analysis
class Constraint_Analysis(Data):
    """Computes the constraint diagram to select initial characteristics

    Assumptions:
    None

    Source:
    None
    """        
    def __defaults__(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """          
        self.tag    = 'Constraint_Analysis'          
        # --------------------------------------------------------------
        #   User inputs
        # --------------------------------------------------------------
        
        self.wingloadinglist_pa = np.arange(5, 125, 5)* Units['force_pound/foot**2'] 
        
        ## Brief 
        # take-off
        self.runway_elevation = 0 #0
        self.ground_run = -1 #1158.24
        # climb
        self.climb_altitude = 0 #0
        self.climb_speed = -1 # 255
        self.climb_rate = -1 #1800         
        # cruise
        self.cruise_altitude = -1 #6096 
        self.cruise_speed = -1 #255
        self.cruise_thrust_fraction = 1.0 #1.0        
        # turn
        self.turn_angle = -1 # 30
        self.turn_load_factor = -1 #1.1547
        self.turn_altitude = 0 #6096
        self.turn_speed = -1 #255
        self.turn_specificenergy = 0 #0
        # ceiling
        self.service_ceiling = -1 #7315.2
        self.service_ceiling_climb_speed = -1 #255          
        # landing
        self.landing_roll = -1 #2100*0.3048
        
        ## Design
        self.aspect_ratio = 8 #14
        self.thickness_to_chord = 0.12 
        #engine
        self.bypass_ratio = -1 #-1
        self.throttle_ratio = 1.07 #1.07
        #sweep        
        self.sweep_leading_edge = 0
        self.sweep_max_thickness = -1
        self.sweep_quarter_chord = -1
        self.oswald_factor = -1
        #high-lift   
        self.high_lift_type_clean = -1 # 0 for clean configuration, we could calculate clmax_clean        
        self.high_lift_type_takeoff = -1
        self.high_lift_type_landing = -1 # 6
        #weight
        self.climb_weight_fraction = 1.0
        self.cruise_weight_fraction = 1.0
        self.turn_weight_fraction = 1.0
        self.ceiling_weight_fraction = 1.0
 
        ## Performance
        # take-off
        self.rolling_resistance = 0.03 #0.05        
        self.cd_takeoff = 0.09 #0.12
        self.cl_takeoff = 0.95 #1.16
        self.cl_max_takeoff = -1 #2.2
        # cruise, turn, ceiling
        self.cd_min_clean = 0.03 #0.0134
        self.cl_max_clean = -1 #1.8999
        # landing
        self.cl_max_approach = -1 #1.8999
        # propeller
        self.propeller_efficiency_takeoff =  0.45 #0.54
        self.propeller_efficiency_climb = 0.75 #0.85
        self.propeller_efficiency_cruise = 0.85 #0.85
        self.propeller_efficiency_turn = 0.85 #0.85
        self.propeller_efficiency_servceil = 0.65 #0.85
        return
    
    
    def create_ADRinput(self):
        """This sets the default values and methods for the analysis.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        designbrief <dictionary>
        designdefinition <dictionary>
        designperformance <dictionary>
        designatm <dictionary>
        
        Properties Used:
        N/A
        """  
        designbrief = {'rwyelevation_m': self.runway_elevation, 
                       'groundrun_m': self.ground_run,
                       'climbalt_m': self.climb_altitude,  
                       'climbspeed_tas': self.climb_speed, 
                       'climbrate_mps': self.climb_rate, 
                       'cruisealt_m': self.cruise_altitude, 
                       'cruisespeed_tas': self.cruise_speed, 
                       'cruisethrustfact': self.cruise_thrust_fraction,
                       'turnangle': self.turn_angle,
                       'stloadfactor': self.turn_load_factor, 
                       'turnalt_m': self.turn_altitude, 
                       'turnspeed_tas': self.turn_speed,  
                       'turnspecificenergy_mps': self.turn_specificenergy, 
                       'servceil_m': self.service_ceiling,
                       'secclimbspd_tas': self.service_ceiling_climb_speed,                       
                       'landingroll_m': self.landing_roll }
                             
        designdefinition = {'aspectratio': self.aspect_ratio,
                            'tc_ratio': self.thickness_to_chord,
                            'bpr': self.bypass_ratio,
                            'tr': self.throttle_ratio,  
                            'sweep_le_rad': self.sweep_leading_edge,  
                            'sweep_mt_rad': self.sweep_max_thickness, 
                            'sweep_025_rad': self.sweep_quarter_chord,
                            'oswaldfactor': self.oswald_factor,
                            'highlift_type':{'clean': self.high_lift_type_clean,
                                             'take-off': self.high_lift_type_takeoff,
                                             'landing': self.high_lift_type_landing},                            
                            'weightfractions':{'climb': self.climb_weight_fraction,
                                               'cruise': self.cruise_weight_fraction,
                                                'turn': self.turn_weight_fraction,                                                
                                                'servceil': self.ceiling_weight_fraction}
                                             } 
                                                        
        designperformance = {'CDTO': self.cd_takeoff, 
                             'CLTO': self.cl_takeoff, 
                             'CLmaxTO': self.cl_max_takeoff, 
                             'mu_R': self.rolling_resistance, 
                             'CDminclean': self.cd_min_clean,
                             'etaprop':{'take-off': self.propeller_efficiency_takeoff, 
                                        'climb': self.propeller_efficiency_climb,
                                        'cruise': self.propeller_efficiency_cruise,
                                        'turn': self.propeller_efficiency_turn, 
                                        'servceil': self.propeller_efficiency_servceil},
                             'CLmaxclean': self.cl_max_clean, 
                             'CLmaxApp': self.cl_max_approach} 
                             
        designatm = at.Atmosphere()
        return designbrief, designdefinition, designperformance, designatm
 
 
    def calculate_thrust_to_weight(self):
        """Calculate the thrust to weight ratio.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        twratio <dictionary>

        Properties Used:

        """         
        concept = ca.AircraftConcept(self.create_ADRinput()[0],self.create_ADRinput()[1], self.create_ADRinput()[2], self.create_ADRinput()[3])          
        thrust_to_weight_ratio = concept.twrequired(self.wingloadinglist_pa,feasibleonly=True)
        landing_wing_loading = concept.wsmaxlanding_pa()
        thrust_to_weight_ratio['landing_ws'] = landing_wing_loading
         
        return thrust_to_weight_ratio
 
 
    def calculate_power_to_weight(self):
        """Calculate the power to weight ratio.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        pwratio <dictionary>

        Properties Used:

        """ 
        concept = ca.AircraftConcept(self.create_ADRinput()[0],self.create_ADRinput()[1], self.create_ADRinput()[2], self.create_ADRinput()[3])        
        power_to_weight_wratio = concept.powertoweightrequired(self.wingloadinglist_pa,feasibleonly=True)
        landing_wing_loading = concept.wsmaxlanding_pa()
        power_to_weight_wratio['landing_ws'] = landing_wing_loading
            
        return power_to_weight_wratio


    def unit_conversion_for_plot(self,thrust_or_power,unit_system):
        """Create the conversions of the units in order to calculate and label the plot.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        thrust_or_power 
        unit_system
        
        Outputs:
        unit_x
        unit_y
        label_x_unit
        label_y_unit
        
        Properties Used:

        """ 
        if thrust_or_power.lower() == 'thrust':
            if unit_system.lower() == "si":
                unit_x = 'newton/meter**2'
                unit_y = 'newton/newton'
                label_x_unit = 'N/m^2'
                label_y_unit = ' - '
            elif unit_system.lower() == "si-kg":
                unit_x = 'force_kilogram /meter**2'
                unit_y = 'newton/newton'
                label_x_unit = 'kg/m^2'
                label_y_unit = ' - '
            elif unit_system.lower() == "us":
                unit_x = 'force_pound/feet**2'
                unit_y = 'force_pound/force_pound'
                label_x_unit = 'lb/ft^2'
                label_y_unit = ' - '
            else:
                error_si_or_us = "Please select 'SI' or 'US' unit system for your constraint diagram (case-insensitive)"
                raise ValueError(error_si_or_us)
            
        elif thrust_or_power.lower() == 'power':
            if unit_system == "si":
                unit_x = 'newton/meter**2'
                unit_y = 'watt/newton'
                label_x_unit = 'N/m^2'
                label_y_unit = 'W/N'
            elif unit_system.lower() == "si-hp/kg":
                unit_x = 'force_kilogram /meter**2'
                unit_y = 'horsepower/force_kilogram'
                label_x_unit = 'kg/m^2'
                label_y_unit = 'hp/kg'
            elif unit_system.lower() == "us":
                unit_x = 'force_pound/feet**2'
                unit_y = 'horsepower/force_pound'
                label_x_unit ='lb/ft^2'
                label_y_unit = 'hp/lbf'
            else:
                error_si_or_us = "Please select 'SI' or 'US' unit system for your constraint diagram (case-insensitive)"
                raise ValueError(error_si_or_us)
        else:
            error_power_or_thrust = "Please select 'Power' or 'Thrust' constraint diagram (case-insensitive)."
            raise ValueError(error_power_or_thrust)
        return unit_x, unit_y, label_x_unit, label_y_unit
        
        
    def plot_constraint_diagram(self, thrust_or_power, unit_system):
        """Plot the constraint diagram and return the design point.

        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:

        """           
 
        if thrust_or_power.lower() == 'thrust':
            constraint_dictionary = self.calculate_thrust_to_weight() 
            label_y = 'Thrust'            
        elif thrust_or_power.lower() == 'power':
            constraint_dictionary = self.calculate_power_to_weight() 
            label_y = 'Power'
        else:
            error_power_or_thrust = "Please select 'Power' or 'Thrust' constraint diagram (case-insensitive)."
            raise ValueError(error_power_or_thrust)
        
        unit_x, unit_y, label_x_unit, label_y_unit = self.unit_conversion_for_plot(thrust_or_power,unit_system)
        
        # unpack constraint dictionary
        take_off = constraint_dictionary['take-off'] / Units[unit_y]
        climb = constraint_dictionary['climb'] / Units[unit_y]
        cruise = constraint_dictionary['cruise'] / Units[unit_y]
        ceiling = constraint_dictionary['servceil'] / Units[unit_y]
        turn = constraint_dictionary['turn'] / Units[unit_y]
        landing = constraint_dictionary['landing_ws'] / Units[unit_x]
        combined = constraint_dictionary['combined'] / Units[unit_y]
        wingloadinglist = self.wingloadinglist_pa / Units[unit_x]     
        
        # determine the design point
        design_point = np.nanmin(combined)
        design_point_ws = wingloadinglist[np.where(combined == np.nanmin(combined))[0][0]]

        if design_point_ws > landing:
            design_point_ws = landing
            index_no_nans = np.logical_not(np.isnan(combined))
            combined_no_nans =  combined[index_no_nans]
            wingloadinglist_no_nans = wingloadinglist[index_no_nans]
            design_point = np.interp(design_point_ws,wingloadinglist_no_nans, combined_no_nans)
       
        print("The calculated design point has {0}-to-Weight ratio {1:f} ({2}) and Wing Loading {3:f} ({4})"\
        .format(label_y, design_point, label_y_unit, design_point_ws, label_x_unit))
        
        
        ##
        # plot the constraint diagram
        plt.rcParams["figure.figsize"] = [10,7]
        plt.plot(wingloadinglist, take_off,  label = 'Take-off')
        plt.plot(wingloadinglist, turn, label = 'Turn')
        plt.plot(wingloadinglist, climb, label = 'Climb')
        plt.plot(wingloadinglist, cruise, label = 'Cruise')
        plt.plot(wingloadinglist, ceiling, label = 'Service ceiling')
        plt.plot(wingloadinglist, combined, label = 'Combined', linewidth=4)
        plt.plot([landing, landing], [0, np.nanmax(combined)], label = 'landing')
        
        plt.plot(design_point_ws, design_point, \
                      label = 'Design Point', marker = 'o', markersize = '10')
            
        legend = plt.legend(loc='upper left')
        plt.ylabel("{0}-to-Weight (${1}$)".format(label_y,label_y_unit))
        plt.xlabel("Wing Loading (${0}$)".format(label_x_unit))
        plt.title("Minimum Sea Level {0}-to-Weight Ratio".format(label_y))
        plt.grid(True)
        plt.show()
        
        return 
     
 
    
