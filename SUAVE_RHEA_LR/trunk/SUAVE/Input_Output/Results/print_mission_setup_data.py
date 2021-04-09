## @ingroup Input_Output-Results
# print_mission_breakdown.py

# Created:  SUAVE team
# Modified: Aug 2020, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
import numpy as np

from SUAVE.Core import Units, Data
import time                     # importing library
import datetime                 # importing library

# ----------------------------------------------------------------------
#  Methods
# ----------------------------------------------------------------------
## @ingroup Input_Output-Results
def print_mission_setup_data(results,filename='mission_setup_data.dat', units="imperial"):
    """This creates a file showing mission information.

    Assumptions:
    None

    Source:
    N/A

    Inputs:
    results.segments.*.conditions.
      frames.
        inertial.position_vector     [m]
        inertial.time                [s]
      aerodynamics.lift_coefficient  [-]
      weights.total                  [kg]
      freestream.  
        mach_number                  [-]
        pressure                     [Pa]
    filename (optional)       <string> Determines the name of the saved file
    units (option)            <string> Determines the type of units used in the output, options are imperial and si

    Outputs:
    filename                  Saved file with name as above

    Properties Used:
    N/A
    """ 

    fid = open(filename,'w')   # Open output file
    fid.write('Output file with the final mission set-up data\n\n') #Start output printing

    for key in results.segments.keys():        #loop for all segments
        
        segment = results.segments[key]

        fid.write(str(segment.tag) + '\n')
        if 'climb' in segment.tag:
            fid.write('Initial altitude:       ' + str(segment.altitude_start) + ' m\n')
            fid.write('Final altitude:         '  + str(segment.altitude_end) + ' m\n')
            if 'air_speed' in segment:
                fid.write('True airspeed:      ' + str(segment.air_speed) + ' m/s \n')
            if 'equivalent_air_speed' in segment:
                fid.write('Equivalent airspeed  ' + str(segment.equivalent_air_speed) + ' m/s \n')
            if 'calibrated_air_speed' in segment:
                fid.write('Calibrated airspeed: ' + str(segment.calibrated_air_speed) + ' m/s \n')
            fid.write('Rate-of-climb: ' + str(segment.climb_rate) + ' m/s \n\n')
        elif 'cruise' in segment.tag:
            fid.write('Altitude: ' + str(segment.altitude) + ' m\n')
            if 'air_speed' in segment:
                fid.write('True airspeed:      ' + str(segment.air_speed) + ' m/s \n')
            if 'equivalent_air_speed' in segment:
                fid.write('Equivalent airspeed  ' + str(segment.equivalent_air_speed) + ' m/s \n')
            if 'calibrated_air_speed' in segment:
                fid.write('Calibrated airspeed: ' + str(segment.calibrated_air_speed) + ' m/s \n')
            if 'mach' in segment:
                fid.write('Mach number: ' + str(segment.mach) + ' \n')
            fid.write('Distance: ' + str(segment.distance/1000) + ' km\n\n')
        elif 'descent' in segment.tag:
            if 'air_speed' in segment:
                fid.write('True airspeed:      ' + str(segment.air_speed) + ' m/s \n')
            if 'equivalent_air_speed' in segment:
                fid.write('Equivalent airspeed  ' + str(segment.equivalent_air_speed) + ' m/s \n')
            if 'calibrated_air_speed' in segment:
                fid.write('Calibrated airspeed: ' + str(segment.calibrated_air_speed) + ' m/s \n')
            if 'mach_number' in segment:
                fid.write('Mach number: ' + str(segment.mach_number) + ' \n') 
            if 'mach_start' in segment:
                 fid.write('Start Mach number: ' + str(segment.mach_start) + ' \n') 
                 fid.write('Finish Mach number: ' + str(segment.mach_end) + ' \n')     
            if 'descent_rate' in segment:
                 fid.write('Rate-of-descent: ' + str(segment.descent_rate) + ' m/s\n') 
            fid.write('End Altitude:' + str(segment.altitude_end) + ' m\n\n')    
            
            


    #done! 
    return

# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')    
