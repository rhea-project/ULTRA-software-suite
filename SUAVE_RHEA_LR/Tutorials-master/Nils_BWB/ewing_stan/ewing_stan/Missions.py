#----------------------------------------------------------------------
#   Imports
# ---------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units
import time


# ----------------------------------------------------------------------
#   Define the Mission
# ----------------------------------------------------------------------
def mission_setup(analyses, vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Mission
    # ------------------------------------------------------------------

    
    mission = SUAVE.Analyses.Mission.Sequential_Segments()
    mission.tag = 'Solar_Flight'

    mission.atmosphere  = SUAVE.Attributes.Atmospheres.Earth.US_Standard_1976()
    mission.planet      = SUAVE.Attributes.Planets.Earth()
    
    # unpack Segments module
    Segments = SUAVE.Analyses.Mission.Segments
    
    # base segment
    base_segment = Segments.Segment()   
    ones_row     = base_segment.state.ones_row
    base_segment.process.iterate.unknowns.network            = vehicle.propulsors.network.unpack_unknowns
    base_segment.process.iterate.residuals.network           = vehicle.propulsors.network.residuals    
    base_segment.process.iterate.initials.initialize_battery = SUAVE.Methods.Missions.Segments.Common.Energy.initialize_battery
    base_segment.state.unknowns.propeller_power_coefficient  = vehicle.propulsors.network.propeller.prop_attributes.Cp  * ones_row(1)/15.
    base_segment.state.residuals.network                     = 0. * ones_row(1)      
    


    
    segment = SUAVE.Analyses.Mission.Segments.Climb.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "climb"

    segment.start_time     = time.strptime("Thu, May 02 10:00:00  2019", "%a, %b %d %H:%M:%S %Y")
    segment.battery_energy = vehicle.propulsors.network.battery.max_energy * 0.95 #Charge the battery to start
    # latitude and longitude of Darmstadt. 
    segment.latitude       = 49.8787   # this defaults to degrees (do not use Units.degrees)
    segment.longitude      = 8.6469 # this defaults to degrees

    segment.state.numerics.number_control_points = 32
   
    # connect vehicle configuration
    segment.analyses.extend(analyses.cruise)  

    ones_row = segment.state.ones_row
    segment.state.unknowns.body_angle = ones_row(1) * 7. * Units.deg    
    segment.altitude_start = 0.0   * Units.km
    segment.altitude_end   = 2.0   * Units.km
    segment.air_speed      = 50.0 * Units['m/s']
    segment.climb_rate     = 4.0   * Units['m/s']

    # add to misison
    mission.append_segment(segment)
    
    
    # ------------------------------------------------------------------    
    #   Cruise Segment: constant speed, constant altitude
    # ------------------------------------------------------------------    
    
    segment = SUAVE.Analyses.Mission.Segments.Cruise.Constant_Speed_Constant_Altitude(base_segment)
    segment.tag = "cruise"
    segment.state.numerics.number_control_points = 32
    
    # connect vehicle configuration
    segment.analyses.extend(analyses.cruise)
    
    # segment attributes     
    segment.altitude       = 2.0  * Units.km 
    segment.air_speed      = 55.0 * Units["m/s"]
    segment.distance       = 500.0 * Units.km
    
    # add to mission
    mission.append_segment(segment)     
    
    
    # ------------------------------------------------------------------
    #   Descent Segment: Constant Speed, Constant Rate
    # ------------------------------------------------------------------

    segment = SUAVE.Analyses.Mission.Segments.Descent.Constant_Speed_Constant_Rate(base_segment)
    segment.tag = "descent"
    segment.state.numerics.number_control_points = 32

    # connect vehicle configuration
    segment.analyses.extend(analyses.cruise)

    ones_row = segment.state.ones_row

    # segment attributes    
    segment.state.unknowns.body_angle = ones_row(1) * 5. * Units.deg
    segment.altitude_end = 0.0   * Units.km
    segment.air_speed    = 55.0  * Units['m/s']
    segment.descent_rate = 0.8   * Units['m/s']

    # add to mission
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
    
    # done!
    return missions 




