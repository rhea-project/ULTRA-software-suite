# Analyses.py
# 
# Created:  Mar 2016, M. Vegh
# Modified: Aug 2017, E. Botero

# ----------------------------------------------------------------------        
#   Imports
# ----------------------------------------------------------------------    

import SUAVE
from SUAVE.Core import Data, Units

# ----------------------------------------------------------------------        
#   Setup Analyses
# ----------------------------------------------------------------------  

def setup(configs):
    
    analyses = SUAVE.Analyses.Analysis.Container()

    # build a base analysis for each config
    for tag,config in configs.items():
        analysis = base(config)
        if tag == 'cruise_spoilers':
            # this is done since drag is not sufficient for the desired profile
            analysis.aerodynamics.settings.spoiler_drag_increment = 0.005       
        analyses[tag] = analysis


    return analyses

# ----------------------------------------------------------------------        
#   Define Base Analysis
# ----------------------------------------------------------------------  

def base(vehicle):

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
    weights = SUAVE.Analyses.Weights.Weights_Tube_Wing()
    weights.vehicle = vehicle
    weights.vehicle.settings = Data()
    weights.vehicle.settings.weight_reduction_factors = Data()
    weights.vehicle.settings.weight_reduction_factors.main_wing = 0.0          # 20% structure weight reduction in main wing 
    weights.vehicle.settings.weight_reduction_factors.fuselage  = 0.0          # 20% structure weight reduction in cabin and aft centerbody
    weights.vehicle.settings.weight_reduction_factors.empennage = 0.0          # 20% structure weight reduction in  HTP, VTP        
    
    analyses.append(weights)

    # ------------------------------------------------------------------
    #  Aerodynamics Analysis
    aerodynamics = SUAVE.Analyses.Aerodynamics.Fidelity_Zero()
    aerodynamics.geometry = vehicle
    aerodynamics.settings.fuselage_xt = 0.0                                 # Fuselage transition location
    aerodynamics.settings.wing_xt     = 0.3                                 # Wing transition location    
    aerodynamics.settings.prop_xt     = 0.3                                 # Propulstion system transition location 
   
    analyses.append(aerodynamics)

    # ------------------------------------------------------------------
    #  Stability Analysis
    stability = SUAVE.Analyses.Stability.Fidelity_Zero()
    stability.geometry = vehicle
    analyses.append(stability)

    # ------------------------------------------------------------------
    #  Energy
    energy= SUAVE.Analyses.Energy.Energy()
    energy.network = vehicle.propulsors
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

    return analyses 
