## @ingroup Input_Output-Results
# print_take_off_results.py 

# Created: S. Karpuk, Oct 2020
# Updated: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
#  no import

# ----------------------------------------------------------------------
#  Print output file with weight breakdown
# ----------------------------------------------------------------------
## @ingroup Input_Output-Results
def print_take_off_results(takeoff,filename = 'takeoff_data.dat'):
    """This creates a file showing take-off results.
    Assumptions:

    Source:

    Inputs:

    Outputs:

    """     
    
    # Imports
    import datetime                 # importing library
    
	# start printing
    fid = open(filename,'w')   # Open output file
    fid.write('Output file with take-off breakdown\n\n') #Start output printing
    fid.write( ' Take-off summary: \n\n')    
    fid.write( ' Take-off field length:          ' + str(takeoff.takeoff_field_length) + ' m\n')
    fid.write( ' Ground roll distance:           ' + str(takeoff.ground_roll) + ' m\n')
    fid.write( ' Maximum lift coefficient:       ' + str(takeoff.maximum_lift_coefficient) + '\n')
    fid.write( ' Take-off lift coefficient:      ' + str(takeoff.lift_coefficient) + '\n')
    fid.write( ' Take-off drag coefficient:      ' + str(takeoff.drag_coefficient) + '\n\n')   

    # Print timestamp
    fid.write('\n'+ 43*'-'+ '\n' + datetime.datetime.now().strftime(" %A, %d. %B %Y %I:%M:%S %p"))
    # done
    fid.close    

# ----------------------------------------------------------------------
#   Module Test
# ----------------------------------------------------------------------
if __name__ == '__main__':
    print(' Error: No test defined ! ')   
