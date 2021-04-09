## @ingroup Methods-Propulsion
# propeller_design.py
# 
# Created:  Mar 2021, S. Karpuk
# Modified: 
#          

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
import numpy as np
from SUAVE.Core  import Units
import matplotlib.pyplot     as plt

# ----------------------------------------------------------------------
#  Propeller Design
# ----------------------------------------------------------------------
    
def propeller_design_LoFi(vehicle,prop_attributes):
    """ Size the propeller based on semi-empirical methods and the cubic spline method.
          
          Inputs:
          Either design power or thrust
          prop_attributes.
            hub radius                       [m]
            tip radius                       [m]
            rotation rate                    [rad/s]
            freestream velocity              [m/s]
            number of blades               
            number of stations
            design lift coefficient
            airfoil data                     

          Outputs:
          Twist distribution                 [array of radians]
          Chord distribution                 [array of meters]
              
          Assumptions:
          1. Cubic spline method
          2. the fub radius is equal to 15% of the disk radius

          Sources:
          1. K. Ibrahim, "Selecting Principal Parameters of Baseline Design Configuration for Twin
             Turboprop Transport Aircraft," in 22nd Applied Aerodynamics Conference and Exhibit,
             Providence, Rhode Island, 2004. AIAA 2004-5069.
          2. Torenbeek, Synthesis of Subsonic AIrcraft Design
          3. Gudmundsson, General Aviation AIrcraft Design

    """   

    # Unpack
    M      = vehicle.cruise_mach
    Mmax   = vehicle.max_cruise_mach

    etac   = prop_attributes.cruise_efficiency
    etah   = prop_attributes.max_cruise_efficiency
    alt    = prop_attributes.design_altitude
    Power  = prop_attributes.design_power       # Design power per propeller
    Nb     = prop_attributes.number_blades
    K_prop = prop_attributes.type               # Twin Low-speed, Twin High-speed, Single Low-speed, Single High-speed

    # Determine the propeller diameter
    if K_prop == 'Twin Low-speed':
        Dp   = 0.232 * (Power/Nb)**0.485
        Mtip = 0.8 
    elif K_prop == 'Twin High-speed':
        Dp = 0.3401 * (Power/Nb)**0.3894
        Mtip = 0.85
    elif K_prop == 'Single Low-speed':
        Dp = 0.776 * (Power/Nb)**0.224
        Mtip = 0.8 
    elif K_prop == 'Single High-speed':
        Dp = 0.984 * (Power/Nb)**0.148
        Mtip = 0.85 

    Dsp = 0.15 * Dp

    # Determine the design propeller RPM
    atmosphere = SUAVE.Analyses.Atmospheric.US_Standard_1976()
    atmo_data = atmosphere.compute_values(alt)

    a   = atmo_data.speed_of_sound[0,0]
    rho = atmo_data.density[0,0]

    n     = a * (Mtip**2-M**2)**0.5/(np.pi*Dp)          # rev/s
    omega = 2*np.pi*n

    # Convert power from kW to W
    PowerW = Power*1e3

    # Determine cubic spline coefficients
    Vc = a * M 
    Vh = a * Mmax

    Asp = 0.25 * np.pi * Dsp**2
    Ap  = 0.25 * np.pi * Dp**2

    Tstatic = 0.85 * PowerW**(2/3)*(2*rho*Ap)**(1/3)*(1-Asp/Ap)

    LHS = np.array([[0,0,0,1],[Vc**3,Vc**2,Vc,1],[3*Vc**2,2*Vc,1,0],[Vh**3,Vh**2,Vh,1]])
    RHS = np.array([Tstatic,etac*PowerW/Vc,-etac*PowerW/Vc**2,etah*PowerW/Vh])

    thrust_spline = np.linalg.solve(LHS, RHS)
    eta_spline    = thrust_spline / PowerW
  
    # Plot the curve
    V_plot   = np.linspace(1, Vh, num=100)
    eta_plot = eta_spline[0]*np.power(V_plot,4)+eta_spline[1]*np.power(V_plot,3)+eta_spline[2]*np.power(V_plot,2)+eta_spline[3]*V_plot 
    #eta_plot = thrust_spline[0]*np.power(V_plot,3)+thrust_spline[1]*np.power(V_plot,2)+thrust_spline[2]*np.power(V_plot,1)+thrust_spline[3] 

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(V_plot , eta_plot, color='tab:blue')
    plt.savefig('eta_p.png')

    # Pack outputs
    prop_attributes.tip_radius             = 0.5*Dp
    prop_attributes.hub_radius             = 0.5*Dsp
    prop_attributes.angular_velocity       = omega
    prop_attributes.RPM                    = omega * 30/np.pi
    prop_attributes.prop_efficiency_coefs  = eta_spline
    

    return prop_attributes
