## @ingroup Methods-Aerodynamics-Xfoil
#run_analysis.py
# 
# Created:  Aug 2019, S. Karpuk
# Modified: 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import sys
import time
import subprocess
import os
from threading                                     import Timer
from SUAVE.Methods.Aerodynamics.XFOIL.purge_files  import purge_files
from SUAVE.Core                                    import redirect

## @ingroup Methods-Aerodynamics-Xfoil
def run_analysis(xfoil_object,results_file,deck_file):
    """ This calls the Xfoil executable and runs an analysis

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        xfoil_object - passed into the  call_avl function  
        
    Outputs:
        results

    Properties Used:
        N/A
    """

    call_xfoil(xfoil_object,deck_file)
    
    return 


def call_xfoil(xfoil_object,in_deck):
    """ This function calls the Xfoil executable and executes analyses

    Assumptions:
        None
        
    Source:
        None

    Inputs:
        xfoil_object

    Outputs:
        exit_status

    Properties Used:
        N/A
    """

    xfoil_regression_flag = xfoil_object.regression_flag
    if xfoil_regression_flag:
        exit_status = 0
    else:
        log_file = xfoil_object.settings.filenames.log_filename
        err_file = xfoil_object.settings.filenames.err_filename
        if isinstance(log_file,str):
            purge_files(log_file)
        if isinstance(err_file,str):
            purge_files(err_file)
        xfoil_call = xfoil_object.settings.filenames.xfoil_bin_name


        ctime = time.ctime() # Current date and time stamp

        with redirect.output(log_file,err_file):
            with open(in_deck,'r') as commands:
                print_output = False

                # Initialize suppression of console window output
                if print_output == False:
                    devnull = open(os.devnull,'w')
                    sys.stdout = devnull       
                    
                # Run AVL
                xfoil_run = subprocess.Popen([xfoil_call],stdout=sys.stdout,stderr=sys.stderr,stdin=subprocess.PIPE)
                timer = Timer(5, xfoil_run.kill)
                for line in commands:
                    xfoil_run.stdin.write(line.encode('utf-8'))
                    xfoil_run.stdin.flush()
                  
                # Terminate suppression of console window output  
                if print_output == False:    
                    sys.stdout = sys.__stdout__                    
                commands.close()
        timer.start()
        exit_status=xfoil_run.communicate()
        
        exit_status = xfoil_run.returncode
        ctime = time.ctime()

    return exit_status

