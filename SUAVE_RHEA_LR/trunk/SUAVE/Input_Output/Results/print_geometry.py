## @ingroup Input_Output-Results
# print_geometry.py

# Modified: Aug 2020, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
import numpy as np

from SUAVE.Core import Units
import time                     # importing library
import datetime                 # importing library

# ----------------------------------------------------------------------
#  Methods
# ----------------------------------------------------------------------
## @ingroup Input_Output-Results
def print_geometry(vehicle,filename='mission_breakdown.dat'):
    """This creates a file showing aircraft geometry

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

    # Output wings data
    n_segments = len(vehicle.wings.keys())
    wing_keys  = list(vehicle.wings.keys()) 

    # write header of file
    fid = open(filename, 'w')  # Open output file
    fid.write('Output file with aircraft geometry data \n\n')
    fid.write('  VEHICLE TAG : ' + vehicle.tag + '\n\n')  
    for i in range(n_segments):
        fid.write(vehicle.wings[wing_keys[i]].tag + ":\n")
        fid.write("AR                     = " + str(vehicle.wings[wing_keys[i]].aspect_ratio) + "\n")
        fid.write("span                   = " + str(vehicle.wings[wing_keys[i]].spans.projected) + " m \n")  
        fid.write("area                   = " + str(vehicle.wings[wing_keys[i]].areas.reference) + " sq m\n") 
        fid.write("root chord             = " + str(vehicle.wings[wing_keys[i]].chords.root) + " m \n")
        fid.write('mean aerodynamic chord ='  + str(vehicle.wings[wing_keys[i]].chords.mean_aerodynamic) + " m \n")
        fid.write("taper ratio            = " + str(vehicle.wings[wing_keys[i]].taper) + "\n")  
        fid.write("quarter-chord sweep    = " + str(vehicle.wings[wing_keys[i]].sweeps.quarter_chord) + " deg \n")   
        fid.write("difedral               = " + str(vehicle.wings[wing_keys[i]].dihedral) + " deg \n")
        fid.write("thickness              = " + str(vehicle.wings[wing_keys[i]].thickness_to_chord) + " \n\n\n") 
        try:
          fid.write("folding penalty           = " + str(vehicle.wings[wing_keys[i]].folding_penalty) + " \n\n\n") 
        except:
          pass     

     # close file
    fid.close  

    return

# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')    
