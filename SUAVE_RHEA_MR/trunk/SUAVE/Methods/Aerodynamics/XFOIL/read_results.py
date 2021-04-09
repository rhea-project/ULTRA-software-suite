## @ingroup Methods-Aerodynamics-Xfoil
#read_results.py
# 
# Created:  Aug 2019, S. Karpuk
# Modified: 


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Units,Data
import os
import numpy as np
## @ingroup Methods-Aerodynamics-Xfoil
def read_results(results_file,flag,Cl,N_panels,N_panel_max):
    """ This functions reads the results from the results text file created 
    at the end of an Xfoil function call

    Assumptions:
        None
        
    Source:

    Inputs:
        results_file              [String]
        flag                      [String]
        Cl                        [Unitless]
        N_panels                  [Unitless]
        N_panel_max               [Unitless]

    Outputs:
        results                   [Data structure with parasitic drag]
        flag                      [String]
        N_panels                  [Unitless]    

    Properties Used:
        N/A
    """    

    results = Data()
    if flag == False:
        try:
            with open(results_file,'r') as res_file:

                lines   = res_file.readlines()
            
                results.alpha_eff = float(lines[12][2:9].strip()) * Units.degrees
                results.Cd_eff    = float(lines[12][19:27].strip())
                results.Cdp_eff   = float(lines[12][28:37].strip())
        except:
            if N_panels != N_panel_max:
                N_panels = N_panels + 20
                flag == False
            else:
                flag = True
    else:
        with open(results_file,'r') as res_file:
            Prelim_results           = Data()
            Prelim_results.alpha_eff = []
            Prelim_results.Cd_eff    = []
            Prelim_results.Cdp_eff   = []
            Prelim_results.Cl        = []
            
            lines = [line.split() for line in res_file]                      

            for i in range(12,len(lines)):
                Prelim_results.alpha_eff.append(float(lines[i][0]) * Units.degrees)
                Prelim_results.Cd_eff.append(float(lines[i][2]))
                Prelim_results.Cdp_eff.append(float(lines[i][3]))
                Prelim_results.Cl.append(float(lines[i][1]))

            Prelim_results.alpha_eff = np.array(Prelim_results.alpha_eff)
            Prelim_results.Cd_eff    = np.array(Prelim_results.Cd_eff )
            Prelim_results.Cdp_eff   = np.array(Prelim_results.Cdp_eff)
            Prelim_results.Cl        = np.array(Prelim_results.Cl)

            if Prelim_results.Cl[len(Prelim_results.Cl)-1] < Cl:
                print(Prelim_results.Cl,Cl)
                flag = True
            else:
                results.alpha_eff = np.interp(Cl,Prelim_results.Cl,Prelim_results.alpha_eff)
                results.Cd_eff    = np.interp(Cl,Prelim_results.Cl,Prelim_results.Cd_eff)
                results.Cdp_eff   = np.interp(Cl,Prelim_results.Cl,Prelim_results.Cdp_eff)
                flag = False


    return results, flag, N_panels
