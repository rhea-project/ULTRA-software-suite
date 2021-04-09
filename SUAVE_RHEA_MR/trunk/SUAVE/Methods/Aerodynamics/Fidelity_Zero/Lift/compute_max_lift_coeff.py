## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
# compute_max_lift_coeff.py
#
# Created:  Dec 2013, A. Variyar
# Modified: Feb 2014, T. Orra
#           Jan 2016, E. Botero        
#           Feb 2019, E. Botero      
#           Jul 2020, E. Botero 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

#SUAVE Imports
import SUAVE
import numpy as np
from SUAVE.Core import Units
from SUAVE.Components import Wings
from SUAVE.Core  import Data

from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Lift.compute_slat_lift import compute_slat_lift
from SUAVE.Methods.Aerodynamics.Fidelity_Zero.Lift.compute_flap_lift import compute_flap_lift

# ----------------------------------------------------------------------
#  compute_max_lift_coeff
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_max_lift_coeff(state,settings,geometry):
    """Computes the maximum lift coefficient associated with an aircraft high lift system

    Assumptions:
    None

    Source:
    Unknown

    Inputs:
    analyses.max_lift_coefficient_factor       [Unitless]
    vehicle.reference_area                     [m^2]
    vehicle.wings.                             
      areas.reference                          [m^2]
      thickness_to_chord                       [Unitless]
      chords.mean_aerodynamic                  [m]
      sweeps.quarter_chord                     [radians]
      taper                                    [Unitless]
      flaps.chord                              [m]
     control_surfaces.flap.deflection          [radians]
     control_surfaces.slat.deflection          [radians]
      areas.affected                           [m^2]
      control_surfaces.flap.configuration_type [string]
    conditions.freestream.                     
      velocity                                 [m/s]
      density                                  [kg/m^3]
      dynamic_viscosity                        [N s/m^2]
                                               
    Outputs:                                   
    Cl_max_ls (maximum CL)                     [Unitless]
    Cd_ind    (induced drag)                   [Unitless]

    Properties Used:
    N/A
    """    


    # initializing Cl and CDi
    Cl_max_ls  = 0
    dfCL0      = 0
    dsCL0      = 0
    vehicle    = geometry
    conditions = state.conditions

    #unpack
    #max_lift_coefficient_factor = settings.maximum_lift_coefficient_factor

    dfCL0   = np.zeros(len(vehicle.wings.keys()))
    alphad0 = np.zeros(len(vehicle.wings.keys()))
    Cla     = np.zeros(len(vehicle.wings.keys()))
    S_flap  = np.zeros(len(vehicle.wings.keys()))
    S_slat  = np.zeros(len(vehicle.wings.keys()))
    c_pr_c  = np.zeros(len(vehicle.wings.keys()))

    i = 0
    for wing in vehicle.wings:

      if wing.vertical == False:
          
        #geometrical data
        Sref       = vehicle.reference_area
        Swing      = wing.areas.reference
        tc         = wing.thickness_to_chord * 100
        chord_mac  = wing.chords.mean_aerodynamic
        sweep_deg  = wing.sweeps.quarter_chord / Units.degree # convert into degrees
        taper      = wing.taper 

        # conditions data
        V    = conditions.freestream.velocity
        roc  = conditions.freestream.density
        nu   = conditions.freestream.dynamic_viscosity

        #--cl max based on airfoil t_c
        Cl_max_ref = 0.9*(-0.0009*tc**3 + 0.0214*tc**2 - 0.053*tc + 0.7005)
        #-reynolds number effect
        Reyn     =  V * roc * chord_mac / nu
        Re_ref   = 9*10**6
        op_Clmax = Cl_max_ref * ( Reyn / Re_ref ) **0.1

        #wing cl_max to outer panel Cl_max
        w_Clmax = op_Clmax* ( 0.919729714285715 -0.044504761904771*taper \
                             -0.001835900000000*sweep_deg +  0.247071428571446*taper**2 +  \
                              0.003191500000000*taper*sweep_deg -0.000056632142857*sweep_deg**2  \
                             -0.279166666666676*taper**3 +  0.002300000000000*taper**2*sweep_deg + \
                              0.000049982142857*taper*sweep_deg**2  -0.000000280000000* sweep_deg**3)

        # Compute CL increment due to Slat
        dsCL0w, dcl_slat, S_slat[i] = compute_slat_lift(wing,w_Clmax)

         # Compute CL increment due to Flap
        dcl_flap, dfCL0w, alphad0[i], c_pr_c[i], S_flap[i], Cla[i] = compute_flap_lift(wing)

        #results
        Cl_max_ls += (w_Clmax + dcl_slat + dcl_flap) * Swing / Sref
        dsCL0     += dsCL0w
        dfCL0[i]   = dfCL0w

      else:
        dfCL0[i]   = 0.0  
        alphad0[i] = 0.0  
      
      i += 1

    Cl_max_ls = Cl_max_ls #* max_lift_coefficient_factor

    return Cl_max_ls, dfCL0, alphad0, c_pr_c, S_flap, S_slat, Cla
