# Procedure.py
# 
# Created:  Apr 2019, N. Kleemann
# Modified: Jul 2019, N. Kleemann

#----------------------------------------------------------------------
#   Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data
#from mpi4py import MPI

import numpy as np


from Vehicles import *
from Analyses import *
from Missions import *
from Plot_Mission import *


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------
def main():

    # build the vehicle, configs, and analyses
    configs, analyses = full_setup()
    
    configs.finalize()
    analyses.finalize() 
    
    # weight analysis
    weights = analyses.configs.base.weights
    breakdown = weights.evaluate()
    print(breakdown)
    analyses.configs.base.mass_properties = weights

    #SUAVE.Vehicle.test_test(analyses.configs.base)
    #mass = SUAVE.Vehicle.sum_mass(analyses.configs.base)
    #print(mass)
    #cg = SUAVE.Vehicle.CG(analyses.configs.base)
    #print(cg)
    
    #for item in analyses.configs.base:
    #    print(type(item))
    #mass = 3.
    #mass, cg = SUAVE.Methods.Center_of_Gravity.compute_mass_and_cg(analyses.configs.base)
    #cg = SUAVE.Methods.Center_of_Gravity.compute_aircraft_center_of_gravity(analyses.configs.base.wings['main_wing'])
    #print(mass, cg)

    #SUAVE.Vehicle.CG(analyses.configs.base)
    
    # stability analysis
    #stability = analyses.configs.base.stability
    #print(stability)

    # mission analysis
    mission = analyses.missions.base
    results = mission.evaluate()

    
    # plot results    
    plot_mission(results)
    
    return

# ----------------------------------------------------------------------
#   Analysis Setup
# ----------------------------------------------------------------------

def full_setup():
    
    # vehicle data
    vehicle  = vehicle_setup()
    configs  = configs_setup(vehicle)
    
    # vehicle analyses
    configs_analyses = analyses_setup(configs)
    
    # mission analyses
    mission  = mission_setup(configs_analyses,vehicle)
    missions_analyses = missions_setup(mission)

    analyses = SUAVE.Analyses.Analysis.Container()
    analyses.configs  = configs_analyses
    analyses.missions = missions_analyses
    
    return configs, analyses


if __name__ == '__main__':
    #comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()
    #name = MPI.Get_processor_name()
    #if rank == 0:
    main()

    plt.show()
