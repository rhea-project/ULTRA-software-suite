## @ingroup Methods-Airspeeds
# Convert_airspeeds.py
# 
# Created:  Aug 2020, S. Karpuk
# Modified: 
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
import numpy as np
from SUAVE.Core import Units
from SUAVE.Analyses.Atmospheric         import US_Standard_1976 as atmosphere

# ----------------------------------------------------------------------
#  Initialize Conditions
# ----------------------------------------------------------------------
## @ingroup Methods-Missions-Segments-Descent
def Convert_airspeeds(index_from,airspeed_in,altitude):
    """Converts a given airspeed into different airspeed types

    Assumptions:

    Source:
    N/A

    Inputs:
    index_from                                  [Unitless]
    airspeed_in                                 [an airspeed value]
    altitude                                    [meters]                                 

    Outputs:
    airspeed_out                                [an airspeed value]

    Properties Used:
    N/A
    """       
    
    # obtain atmospheric data
    atmo      = atmosphere()
    atmo_data = atmo.compute_values(altitude)
    p  = atmo_data.pressure[0,0]
    T  = atmo_data.temperature[0,0]
    ro = atmo_data.density[0,0]
    a  = atmo_data.speed_of_sound[0,0]
    mu = atmo_data.dynamic_viscosity[0,0]

    atmo_SL      = atmosphere()
    atmo_data_SL = atmo_SL.compute_values(0)
    p0  = atmo_data_SL.pressure[0,0]
    ro0 = atmo_data_SL.density[0,0]



    if (index_from == "tas" or index_from == "TAS"):
        # Convert from TAS
        tas  = airspeed_in
        eas  = tas * np.sqrt(ro/ro0)
        mach = tas / a 

        qc  = p * ((1+0.2*mach)**3.5-1)
        cas = eas * np.sqrt(p0/p) * np.sqrt(((qc/p0+1)**0.286-1)/((qc/p+1)**0.286-1)) 

    if (index_from == "cas" or index_from == "CAS"):
        # Convert from CAS
        cas = airspeed_in
        cas_US = cas / Units.knots
        mach   = 2.236*np.sqrt((((1+4.575*10**(-7)*cas_US**2)**3.5-1)/(p/p0)+1)**0.2857-1)
        qc     = p * ((1+0.2*mach)**3.5-1)
        eas    = cas * np.sqrt(p/p0) * np.sqrt(((qc/p+1)**0.286-1)/((qc/p0+1)**0.286-1)) 
        tas    = eas/np.sqrt(ro/ro0)

    if (index_from == "eas" or index_from == "EAS"):
        # Convert from EAS
        eas = airspeed_in
        tas = eas/np.sqrt(ro/ro0)
        mach = tas / a 
        qc   = p * ((1+0.2*mach)**3.5-1)       
        cas  = eas * np.sqrt(p0/p) * np.sqrt(((qc/p0+1)**0.286-1)/((qc/p+1)**0.286-1)) 

    if (index_from == "mach" or index_from == "Mach"):
        # Convert from Mach
        mach = airspeed_in
        tas  = mach * a
        eas  = tas * np.sqrt(ro/ro0) 
        qc   = p * ((1+0.2*mach)**3.5-1)
        cas  = eas * np.sqrt(p0/p) * np.sqrt(((qc/p0+1)**0.286-1)/((qc/p+1)**0.286-1)) 

    # Pack outputs
    airspeeds = np.array([cas, eas, tas, mach])


    return airspeeds

