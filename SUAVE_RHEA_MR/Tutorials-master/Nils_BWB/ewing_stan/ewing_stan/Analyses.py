# Procedure.py
# 
# Created:  Apr 2019, N. Kleemann
# Modified: July 2019, N. Kleemann


# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import SUAVE


# ----------------------------------------------------------------------
#   Define the Vehicle Analyses
# ----------------------------------------------------------------------

def analyses_setup(configs):
    
    analyses = SUAVE.Analyses.Analysis.Container()
    
    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base_analysis(config)
        analyses[tag] = analysis
    
    return analyses

def base_analysis(vehicle):

    # ------------------------------------------------------------------
    #   Initialize the Analyses
    # ------------------------------------------------------------------     
    analyses = SUAVE.Analyses.Vehicle()
    
    # ------------------------------------------------------------------
    #  Basic Geometry Relations
    sizing = SUAVE.Analyses.Sizing.Sizing()
    sizing.features.vehicle = vehicle
    analyses.append(sizing)
    
    # ------------------------------------------------------------------
    #  Weights
    weights = SUAVE.Analyses.Weights.Weights_UAV()
    #weights = SUAVE.Analyses.Weights.Weights_Electric_Stopped_Rotor()
    #weights = SUAVE.Analyses.Weights.Weights_Flying_Wing()
    weights.vehicle = vehicle
    analyses.append(weights)
    
    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.AVLQ3D()
    # aerodynamics = SUAVE.Analyses.Aerodynamics.AVL()
    aerodynamics.settings.spanwise_vortices = 15
    aerodynamics.settings.chordwise_vortices = 10
    aerodynamics.geometry = vehicle
    aerodynamics.settings.drag_coefficient_increment = 0.000
    analyses.append(aerodynamics)
    
    # ------------------------------------------------------------------
    #  Energy
    energy = SUAVE.Analyses.Energy.Energy()
    energy.network = vehicle.propulsors #what is called throughout the mission (at every time step))
    analyses.append(energy)
    
    # ------------------------------------------------------------------
    #  Planet Analysis
    planet = SUAVE.Analyses.Planets.Planet()
    analyses.append(planet)
    
    # ------------------------------------------------------------------
    #  Atmosphere Analysis
    atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmosphere.features.planet = planet.features
    analyses.append(atmosphere) 

    # ------------------------------------------------------------------
    #  Stability Analysis
    #stability = SUAVE.Analyses.Stability.AVL()
    #stability.geometry = vehicle
    #analyses.append(stability) 

    
    # done!
    return analyses  
