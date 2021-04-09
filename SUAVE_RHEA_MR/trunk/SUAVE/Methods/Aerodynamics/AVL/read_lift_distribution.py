## @ingroup methods-Aerodynamics-AVL
# read_lift_distribution.py
#
# Created: Jul 2019, S. Karpuk 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data

import glob
import re
import numpy as np


def read_lift_distribution(AoA,Mach,Y_strip):
    # Search for Trefftz plane files in the directory
    # and read the data

    Lift_distribution = SUAVE.Core.ContainerOrdered()

    files = glob.glob('results*')
    files.sort(key = alphanum_key)
    
    # Attach a string exstention to the AoA/Mach case     
    for ii in range(len(files)):
        with open(files[ii],'r') as res_file:
            case        = Data()
            case.Y      = []
            case.Ystrip = []
            case.chord  = []
            case.area   = []
            case.Cl     = []
            case.ClQ3D  = []
            lines       = res_file.readlines()
            case.tag    = lines[13][11:25].strip()
            case.AoA    = float(lines[15][11:20].strip())
            case.Mach   = float(lines[17][11:20].strip())

        # Read the Trefftz plane file
        for results_filename in glob.glob('Trefftz_results_' + case.tag + '*'):
            with open(results_filename,'r') as res_file1:
                lines = res_file1.readlines()

            # Search for the wing lift distribution and obtain the data
            with open(results_filename,'r') as res_file1:
                line_index = 0
                for line  in res_file1:
                    for part in line.split():        
                        if "main_wing" in part:
                            i = 1
                            while lines[line_index+13+i][8:10].strip() != "" and lines[line_index+13+i][8:10].strip() != "--":
                                case.Y.append(float(lines[line_index+13+i][7:16].strip()))                           
                                case.chord.append(float(lines[line_index+13+i][16:25].strip()))
                                case.area.append(float(lines[line_index+13+i][26:34].strip()))
                                case.Cl.append(float(lines[line_index+13+i][64:70].strip()))
                                i = i + 1
                    line_index = line_index + 1
        case.Y      = np.array(case.Y)
        case.chord  = np.array(case.chord)
        case.area   = np.array(case.area)
        case.Cl     = np.array(case.Cl)
        
        argsort     = case.Y.argsort()
        case.Y      = case.Y[argsort]
        case.chord  = case.chord[argsort]
        case.area   = case.area[argsort]
        case.Cl     = case.Cl[argsort]
        
        # Interpolate the lift distribution for the Q3D
        case.Ystrip =Y_strip 
        case.ClQ3D  = np.interp(Y_strip,case.Y,case.Cl)
        Lift_distribution.append(case)   

    return Lift_distribution
                    

def tryint(s):
    try:  
        return int(s)
    except:
        return s

def alphanum_key(s):
    return[tryint(c) for c in re.split('([0-9]+)',s)]

    
    
