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
        fid.write('mean aerodynamic chord = ' + str(vehicle.wings[wing_keys[i]].chords.mean_aerodynamic) + " m \n")
        fid.write("taper ratio            = " + str(vehicle.wings[wing_keys[i]].taper) + "\n")  
        fid.write("quarter-chord sweep    = " + str(vehicle.wings[wing_keys[i]].sweeps.quarter_chord/Units.degrees) + " deg \n")  
        fid.write('leading edge sweep     = ' + str(vehicle.wings[wing_keys[i]].sweeps.leading_edge/Units.degree)   + ' deg \n') 
        fid.write("dihedral               = " + str(vehicle.wings[wing_keys[i]].dihedral/Units.degrees) + " deg \n")
        fid.write("thickness              = " + str(vehicle.wings[wing_keys[i]].thickness_to_chord) + "\n") 
        try:
          fid.write("folding penalty        = " + str(vehicle.wings[wing_keys[i]].folding_penalty) + " \n\n") 
        except:
          fid.write("\n\n") 
 
        try: 
          fid.write("Tail volume ratio      = " + str(vehicle.wings[wing_keys[i]].volume_ratio) + " \n\n\n") 
        except:
          fid.write("\n\n")    

        #try:
        n_sections = vehicle.wings[wing_keys[i]].Segments.keys()
        for j in range(len(n_sections)):
          fid.write(str(vehicle.wings[wing_keys[i]].Segments[j].tag) + "\n")    
          fid.write("percent span location               = " + str(vehicle.wings[wing_keys[i]].Segments[j].percent_span_location) + "\n") 
          fid.write("twist                               = " + str(vehicle.wings[wing_keys[i]].Segments[j].twist/Units.degrees)   + " deg \n") 
          fid.write("root chord percent                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].root_chord_percent)    + " \n") 
          fid.write("dihedral                            = " + str(vehicle.wings[wing_keys[i]].Segments[j].dihedral_outboard/Units.degree)     + " deg \n") 
          fid.write("quarter-chord sweep                 = " + str(vehicle.wings[wing_keys[i]].Segments[j].sweeps.quarter_chord/Units.degree)  + " deg \n") 
          fid.write("thickness_to_chord                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].thickness_to_chord )  + " \n") 
          fid.write("front spar percent                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].front_spar )  + " \n") 
          fid.write("rear spar percent                   = " + str(vehicle.wings[wing_keys[i]].Segments[j].rear_spar )  + " \n") 
          fid.write("transition_x_upper                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].transition_x_upper )  + " \n") 
          fid.write("transition_x_lower                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].transition_x_lower )  + " \n") 
          fid.write("transition_x_lower                  = " + str(vehicle.wings[wing_keys[i]].Segments[j].transition_x_lower )  + " \n\n") 

        #except:
        fid.write("\n\n")   


     # close file
    fid.close  

    return

# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')    
