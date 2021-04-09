## @ingroup Analyses-Constraint_Analysis-ADR

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Modified: Aug 2020, R.Y. Yanev 

"""
.. _constraints_module:

Constraint Analysis Module
--------------------------

This module contains tools for the constraint analysis of fixed
wing aircraft.

"""


__author__ = "Andras Sobester"

# pylint: disable=locally-disabled, too-many-instance-attributes
# pylint: disable=locally-disabled, too-many-branches
# pylint: disable=locally-disabled, too-many-statements
# pylint: disable=locally-disabled, too-many-locals
# pylint: disable=locally-disabled, too-many-lines

import math
import warnings
from scipy import constants
import numpy as np

from SUAVE.Analyses.Constraint_Analysis.ADRpy import atmospheres as at
from SUAVE.Analyses.Constraint_Analysis.ADRpy import unitconversions as co
from SUAVE.Analyses.Constraint_Analysis.ADRpy import mtools4acdc as actools


class AircraftConcept:
    """Definition of a basic aircraft concept. **This is the most important
    class in ADRpy.** An object of this class defines an
    aircraft design in terms of the *brief* it is aiming to meet, high
    level *design* variables that specify it, key parameters that describe
    its *performance*, as well as the *atmosphere* it operates in. These are
    the four arguments that define an object of the AircraftConcept class.
    The first three are dictionaries, as described below, the last is an object
    of `Atmosphere <https://adrpy.readthedocs.io/en/latest/#atmospheres.Atmosphere>`_
    class.

    **Parameters:**

    brief
        Dictionary. Definition of the design brief, that is, the requirements
        the design seeks to meet. Contains the following key names:

        rwyelevation_m
            Float. The elevation (in metres) of the runway againts which the take-off
            constraint is defined. Optional, defaults to zero (sea level).

        groundrun_m
            Float. Length (in metres) of take-off ground run in meters at the elevation
            defined by the *rwyelevation_m* entry of the dictionary. This is a basic,
            100% N1, no wind, zero runway gradient ground run.
     
        landingroll_m
            Float. Lenght (in metres) of landing roll distance.   
        
        stloadfactor
            Float. Load factor to be sustained by the aircraft in a steady, level turn.

        turnalt_m
            Float. Altitude (in metres) where the turn requirement is defined.
            Optional, defaults to zero (sea level).

        turnspeed_tas
            Float. True airspeed (in meters per second) at which the turn requirement (above) has to be met.
            Since the dynamics of turning flight is dominated by inertia, which depends
            on ground speed, the turn speed is specified here as TAS (on the zero wind assumption).
            If you'd rather specify this as IAS/CAS/EAS,
            use `eas2tas <https://adrpy.readthedocs.io/en/latest/#atmospheres.Atmosphere.eas2tas>`_
            first to obtain the TAS value.
            
        turnspecificenergy_mps
            Float. The turn specific energy density (in metres per second).
            Optional, defaults to zero (sustained turn).
        
        climbalt_m
            Float. The altitude (in metres) where the climb rate requirement is specified.
            Optional, defaults to zero (sea level).

        climbspeed_tas
            Float. The airspeed (in meters per second, true) at which the required climb
            rate has to be achieved.

        climbrate_mps
            Float. Required climb rate (in meters per second) at the altitude specified
            in the *climbalt_m* entry (above).

        cruisealt_m
            Float. The altitude at which the cruise speed requirement will be defined.

        cruisespeed_tas
            Float. The required cruise speed (in meters per second, true airspeed) at the
            altitude specified in the *cruisealt_m* entry (above).

        cruisethrustfact
            Float. The fraction (nondimensional) of the maximum available thrust at which
            the cruise speed requirement must be achieved.

        servceil_m
            Float. The required service ceiling in meters (that is, the altitude at which
            the maximum rate of climb drops to 100 feet per minute).

        secclimbspd_tas
            Float. The speed (knots true airspeed) at which the service ceiling
            must be reached. This should be an estimate of the best rate of climb speed.

    design
        Dictionary. Definition of key, high level design variables that define the future
        design.

        aspectratio
            Float. Wing aspect ratio.
            
        tc_ratio
            Float. Wing thickness-to-chord ratio.
            
        sweep_le_rad 
            Float. Main wing leading edge sweep angle (in radians). Optional, defaults to
            zero (no sweep).

        sweep_mt_rad
            Float. Main wing sweep angle measured at the maximum thickness point (in radians). Optional,
            defaults to zero.
                   
        sweep_025_rad
            Float. Main wing sweep angle measured at the 25 percent of the chord (in radians). Optional,
            defaults to zero.
     
        highlift_type
            Dictionary, specifying the type of high-lift device for take-off and landing. It 
            should contain the following keys: *take-off* and *landing*. The values for each of the 
            keys can be set to 0 for no flaps, 1 for plain flaps, 2 for single-slotted flaps, 
            3 for Fowler flaps, 4 for double-slotted flaps, 5 for double slotted Fowler flaps
            with slats, 6 for tripple-slotted Fowler flaps with slats.
      
        bpr
            Float. Specifies the propulsion system type. For jet engines (powered by axial
            gas turbines) this should be the bypass ratio (hence *'bpr'*). Set to -1 for
            piston engines, -2 for turboprops and -3 if no power/thrust corrections are needed
            (e.g., for electric motors).

        tr
            Float. Throttle ratio for gas turbine engines. *tr = 1* means that the Turbine Entry
            Temperature will reach its maximum allowable value in sea level standard day
            conditions, so higher ambient temperatures will result in power loss. Higher *tr*
            values mean thrust decay starting at higher altitudes.

        weightfractions
            Dictionary, specifying at what fraction of the maximum take-off weight do various
            constraints have to be met. It should contain the following keys: *take-off*,
            *climb*, *cruise*, *turn*, *servceil*. Optional, each defaults to 1.0 if not
            specified.

    performance
        Dictionary. Definition of key, high level design performance estimates.

        CDTO
            Float. Take-off drag coefficient.

        CLTO
            Float. Take-off lift coefficient.

        CLmaxTO
            Float. Maximum lift coefficient in take-off conditions.

        CLmaxclean
            Float. Maximum lift coefficient in flight, in clean configuration.
        
        CLmaxApp
            Float. Maximum lift coefficient in landing configuration.
       
        mu_R
            Float. Coefficient of rolling resistance on the wheels.

        CDminclean
            Float. Zero lift drag coefficient in clean configuration.

        etaprop
            Dictionary. Propeller efficiency in various phases of the mission.
            It should contain the following keys: *take-off*, *climb*, *cruise*,
            *turn*, *servceil*. Optional, unspecified entries in the dictionary
            default to the following values::

                etap = {'take-off': 0.45, 'climb': 0.75, 'cruise': 0.85,
                        'turn': 0.85, 'servceil': 0.65}

    designatm
            `Atmosphere <https://adrpy.readthedocs.io/en/latest/#atmospheres.Atmosphere>`_
            class object. Specifies the virtual atmosphere in which all the design
            calculations within the *AircraftConcept* class will be performed.
    """

    def __init__(self, brief, design, performance, designatm):

        ## Assign a default, if needed, to the atmosphere
        if not designatm:
            designatm = at.Atmosphere()
        self.designatm = designatm    
        
        
        ## Unpack the design brief dictionary first:
        # take-off
        self.groundrun_m = brief['groundrun_m']       
        self.rwyelevation_m = brief['rwyelevation_m']  
        
        # climb
        self.climbalt_m = brief['climbalt_m']
        self.climbrate_mps = brief['climbrate_mps']       
        if design['bpr'] == -1:
            self.climbspeed_tas = 9999
        else: 
            self.climbspeed_tas = brief['climbspeed_tas']     
            
        # cruise
        self.cruisealt_m = brief['cruisealt_m']        
        self.cruisespeed_tas = brief['cruisespeed_tas']
        self.cruisethrustfact = brief['cruisethrustfact']    
        
        # turn
        self.turnalt_m = brief['turnalt_m']
        self.turnspeed_tas = brief['turnspeed_tas']     
        self.turnspecificenergy_mps = brief['turnspecificenergy_mps']     
        if brief['stloadfactor'] == -1 and brief['turnangle'] !=-1:
            self.stloadfactor = 1/math.cos(brief['turnangle'])
        else:
            self.stloadfactor = brief['stloadfactor']       
            
        # ceiling
        self.servceil_m = brief['servceil_m']
        if design['bpr'] == -1:
            self.secclimbspd_tas = 9999
        else: 
            self.secclimbspd_tas = brief['secclimbspd_tas'] 
            
        # landing
        self.landingroll_m = brief['landingroll_m']
        
        
        ## Unpack the design dictionary next:
        self.aspectratio = design['aspectratio']
        self.oswaldfactor = design['oswaldfactor']
        self.tc_ratio = design['tc_ratio']
        # engine
        self.bpr = design['bpr']
        self.throttle_r = design['tr']
        
        # sweep
        self.sweep_le_rad = design['sweep_le_rad']
        self.sweep_le_deg = math.degrees(self.sweep_le_rad)
        if design['sweep_mt_rad'] != -1:
            self.sweep_mt_rad = design['sweep_mt_rad']
        else:
            self.sweep_mt_rad = self.sweep_le_rad  
            
        if design['sweep_025_rad'] != -1:
            self.sweep_025_rad = design['sweep_025_rad']
            self.sweep_025_deg = math.degrees(self.sweep_025_rad)
        else: 
            self.sweep_025_rad = self.sweep_le_rad
            self.sweep_025_deg = math.degrees(self.sweep_025_rad)
        
        # high-lift
        if design['highlift_type']['take-off'] in np.arange(0,7):
            self.highlift_takeoff = design['highlift_type']['take-off']
        else:
            # Flag if not specified, error thrown by take-off constraint 
            # if CLmaxTO also not specified 
            self.highlift_takeoff = -1   

        if design['highlift_type']['landing'] in np.arange(0,7):  
            self.highlift_landing = design['highlift_type']['landing']           
        else:
            # Flag if not specified, error thrown by landing constraint 
            # if CLmaxApp also not specified
            self.highlift_landing = -1
                    
        if design['highlift_type']['clean'] in np.arange(0,7):
            self.highlift_clean = design['highlift_type']['clean']
        else:
            # Flag if not specified, error thrown by landing constraint 
            # if CLmaxApp also not specified
            self.highlift_clean = -1        
              
        # weight
        self.climb_weight_fraction = design['weightfractions']['climb']
        self.cruise_weight_fraction = design['weightfractions']['cruise']
        self.turn_weight_fraction = design['weightfractions']['turn']
        self.sec_weight_fraction = design['weightfractions']['servceil']
                

        ## Next, unpack the performance dictionary
        # take-off
        self.mu_r = performance['mu_R']
        self.cdto = performance['CDTO']
        self.clto = performance['CLTO']
        self.clmaxto = performance['CLmaxTO']
        
        # cruise, turn, ceiling
        self.cdminclean = performance['CDminclean']
        self.clmaxclean = performance['CLmaxclean']
        
        # landing
        self.clmaxapproach = performance['CLmaxApp']
        
        # propeller
        self.etaprop_to = performance['etaprop']['take-off']
        self.etaprop_climb = performance['etaprop']['climb']
        self.etaprop_cruise = performance['etaprop']['cruise']
        self.etaprop_turn = performance['etaprop']['turn']
        self.etaprop_sec = performance['etaprop']['servceil']       

    # Three different estimates the Oswald efficieny factor:

    def oswaldspaneff1(self):
        """Raymer's Oswald span efficiency estimate, straight wings, moderate AR"""
        return 1.78 * (1 - 0.045 * (self.aspectratio ** 0.68)) - 0.64

    def oswaldspaneff2(self):
        """Oswald span efficiency estimate due to Brandt et al."""
        sqrtterm = 4 + self.aspectratio ** 2 * (1 + (math.tan(self.sweep_mt_rad)) ** 2)
        return 2/(2 - self.aspectratio + math.sqrt(sqrtterm))

    def oswaldspaneff3(self):
        """Raymer's Oswald span efficiency estimate, swept wings > 30 degrees"""
        return 4.61 * (1 - 0.045 * (self.aspectratio ** 0.68)) * \
        ((math.cos(self.sweep_le_rad)) ** 0.15) - 3.1

    def raymerspanneffsweep0to30(self):
        """Raymer's Oswald span efficiency estimate between 0 and 30 sweep"""
        return (1-self.sweep_le_deg/30)*self.oswaldspaneff1() + self.sweep_le_deg/30*self.oswaldspaneff3()

    def induceddragfact(self):
        """Lift induced drag factor k estimate (Cd = Cd0 + k.Cl^2)"""       
        # k = 1 / pi.AR.e       
        if self.oswaldfactor == -1 :
            if abs(self.sweep_le_deg) >= 0 and abs(self.sweep_le_deg) <= 30:
                oswaldspaneff = 0.5 * (self.oswaldspaneff2() + self.raymerspanneffsweep0to30())
            else:
                oswaldspaneff = 0.5 * (self.oswaldspaneff2() + self.oswaldspaneff3())
        else:
            oswaldspaneff = self.oswaldfactor
            
        return 1.0 / (math.pi * self.aspectratio * oswaldspaneff)


    def bestclimbspeedprop(self, wingloading_pa, altitude_m):
        """The best rate of climb speed for a propeller aircraft"""

        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)
        dragfactor = np.sqrt(self.induceddragfact() / (3 * self.cdminclean))
        densfactor = 2 / self.designatm.airdens_kgpm3(altitude_m)

        # Gudmundsson, eq. (18-27)
        bestspeed_mps = np.sqrt(densfactor * wingloading_pa * dragfactor)

        if len(bestspeed_mps) == 1:
            return bestspeed_mps[0]

        return bestspeed_mps


    def thrusttoweight_takeoff(self, wingloading_pa):
        """The thrust to weight ratio required for take-off. This function is an
        implementation of the following simple, analytical model:

        .. math::

            \\frac{\\overline{T}}{W} = 1.21\\frac{W/S}{\\rho C_\\mathrm{Lmax}^\\mathrm{TO}gd_
            \\mathrm{G}}+\\frac{1}{2}\\frac{C_\\mathrm{D}^\\mathrm{TO}}{C_\\mathrm{L}^\\mathrm{TO}}
            +\\frac{1}{2}\\mu_\\mathrm{R}

        where :math:`\\overline{T}` is the average thrust during the take-off run,
        :math:`W/S` is the wing loading, :math:`d_\\mathrm{G}` is the required ground
        roll, :math:`C_\\mathrm{D}^\\mathrm{TO}` and :math:`C_\\mathrm{L}^\\mathrm{TO}`
        are the 'all wheels on the runway' drag and lift coefficient respectively
        in the take-off configuration, :math:`C_\\mathrm{Lmax}^\\mathrm{TO}` is the maximum
        lift coefficient achieved during the take-off run (during rotation), :math:`\\rho`
        is the ambient density and :math:`\\mu_\\mathrm{R}` is the coefficient of rolling
        resistance on the wheels.

        This is a function exposed to the user for clarity and added flexibility.
        If you need to calculate the thrust to weight ratio required for take-off, use
        ``twrequired_to``. This corrects the output of this function to account for the
        environmental conditions (including their impact on engine performance) and includes
        a mapping to static thrust. ``thrusttoweight_takeoff`` should only be used if you
        would like to perform these corrections in a different way than implemented in
        ``twrequired_to``.

        If a full constraint analysis is required, ``twrequired`` should be used.
        A similar 'full constraint set' function is available for calculating the
        power demanded of the engine or electric motor of a propeller-driven aircraft
        (to satisfy the constraint set) - this is called ``powerrequired``.
        """

        groundrun_m = self.groundrun_m

        # Assuming that the lift-off speed is equal to VR, which we estimate at 1.1VS1(T/O)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.rwyelevation_m)          
        
        if self.clmaxto == -1:
            if self.highlift_takeoff != -1:
                self.clmaxto = self.clmax_highlift(self.highlift_takeoff)
            else:
                clmaxmsg = "CLmaxTO must be specified in the performance dictionary or highlift_type 'take-off' in the design dictionary."
                raise ValueError(clmaxmsg)         
        
        vs1to_mps = np.sqrt((2 * wingloading_pa) / (density_kgpm3 * self.clmaxto))

        liftoffspeed_mpstas = 1.1 * vs1to_mps       
        q_liftoff = 1.21 * wingloading_pa/(2 * self.clmaxto)    
        
        thrusttoweightreqd = (liftoffspeed_mpstas ** 2) / (2 * constants.g * groundrun_m) + \
        q_liftoff * self.cdto / wingloading_pa + \
        self.mu_r*(1-q_liftoff*self.clto/wingloading_pa)
        
        return thrusttoweightreqd, liftoffspeed_mpstas


    def thrusttoweight_turn(self, wingloading_pa):
        """Thrust to weight ratio required for a turn at certain altitude, given load factor 
        and specific energy density."""

        nturn = self.stloadfactor
        turnalt_m = self.turnalt_m
        turnspeed_mps = self.turnspeed_tas
        turnspecificenergy_mps = self.turnspecificenergy_mps

        qturn = self.designatm.dynamicpressure_pa(airspeed_mps=turnspeed_mps, altitudes_m=turnalt_m)        
        inddragfact = self.induceddragfact()
        cdmin = self.cdminclean

        # calculate compressibility drag
        mach = self.designatm.mach(turnspeed_mps, turnalt_m)       
        clrequired = nturn * wingloading_pa / qturn        
        cd_c = self.compressibility_drag_wing(mach, clrequired, self.tc_ratio, self.sweep_025_rad)
        
        # calculate twratio       
        twreqtrn = qturn * ((cdmin + cd_c) / wingloading_pa + inddragfact * ((nturn / qturn) ** 2) * wingloading_pa) +\
        turnspecificenergy_mps/turnspeed_mps
        
        return twreqtrn, clrequired


    def _altcorr(self, temp_c, pressure_pa, mach, density_kgpm3):
        """Altitude corrections, depending on propulsion system type"""
        if self.bpr == -1:
            twratio_altcorr = at.pistonpowerfactor(density_kgpm3)
        elif self.bpr == -2:
            twratio_altcorr = at.turbopropthrustfactor(temp_c, pressure_pa, mach, \
            self.throttle_r)
        elif self.bpr == -3: # no correction required
            twratio_altcorr = 1
        elif self.bpr == 0:
            twratio_altcorr = at.turbojetthrustfactor(temp_c, pressure_pa, mach, \
            self.throttle_r, False)
        elif self.bpr < 5:
            twratio_altcorr = at.turbofanthrustfactor(temp_c, pressure_pa, mach, \
            self.throttle_r, "lowbpr")
        else:
            twratio_altcorr = at.turbofanthrustfactor(temp_c, pressure_pa, mach, \
            self.throttle_r, "highbpr")
        return twratio_altcorr


    def twrequired_to(self, wingloading_pa):
        """Calculate the T/W required for take-off for a range of wing loadings

        **Parameters**

        wingloading_pa
            float or numpy array, list of wing loading values in Pa.

        **Returns**

        twratio
            array, thrust to weight ratio required for the given wing loadings.

        liftoffspeed_mpstas
            array, liftoff speeds (TAS - true airspeed) in m/s.

        avspeed_mpstas
            average speed (TAS) during the take-off run, in m/s.

        **See also** ``twrequired``

        **Notes**

        1. The calculations here assume a 'no wind' take-off, conflating ground speed (GS) and
        true airspeed (TAS).

        2. Use `twrequired` if a full constraint analysis is desired, as this integrates
        the take-off, turn, climb, cruise, and service ceiling constraints, as well as
        computing the combined constraint boundary.
        """
        if self.groundrun_m == -1:
            tomsg = "Ground run not specified in the designbrief dictionary."
            raise ValueError(tomsg)

        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)

        twratio, liftoffspeed_mpstas = self.thrusttoweight_takeoff(wingloading_pa)

        # What does this required T/W mean in terms of static T/W required?
        twratio = self.map2static() * twratio

        # What SL T/W will yield the required T/W at the actual altitude?
        temp_c = self.designatm.airtemp_c(self.rwyelevation_m)
        pressure_pa = self.designatm.airpress_pa(self.rwyelevation_m)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.rwyelevation_m)

        for i, los_mps in enumerate(liftoffspeed_mpstas):
            mach = self.designatm.mach(los_mps, self.rwyelevation_m)
            corr = self._altcorr(temp_c, pressure_pa, mach, density_kgpm3)
            twratio[i] = twratio[i] / corr
            
        avspeed_mpstas = liftoffspeed_mpstas / np.sqrt(2)

        if len(twratio) == 1:
            return twratio[0], liftoffspeed_mpstas[0], avspeed_mpstas[0]

        return twratio, liftoffspeed_mpstas, avspeed_mpstas


    def twrequired_trn(self, wingloading_pa):
        """Calculates the T/W required for turning for a range of wing loadings

        **Parameters**

        wingloading_pa
            float or numpy array, list of wing loading values in Pa.

        **Returns**

        twratio
            array, thrust to weight ratio required for the given wing loadings.

        clrequired
            array, lift coefficient values required for the turn (see notes).

        feasibletw
            as twratio, but contains NaNs in lieu of unachievable (CLmax exceeded) values.

        **See also** ``twrequired``

        **Notes**

        1. Use `twrequired` if a full constraint analysis is desired, as this integrates
        the take-off, turn, climb, cruise, and service ceiling constraints, as well as
        computing the combined constraint boundary.

        2. At the higher end of the wing loading range (low wing area values) the CL required
        to achieve the required turn rate may exceed the maximum clean CL (as specified in the
        `CLmaxclean` entry in the `performance` dictionary argument of the `AircraftConcept`
        class object being used). This means that, whatever the T/W ratio, the wings will stall
        at this point. The basic T/W value will still be returned in `twratio`, but there is
        another output, `feasibletw`, which is an array of the same T/W values, with those
        values blanked out (replaced with NaN) that cannot be achieved due to CL exceeding
        the maximum clean lift coefficient.
        """

        if self.turnspeed_tas == -1:
            turnmsg = "Turn speed not specified in the designbrief dictionary."
            raise ValueError(turnmsg)

        if self.stloadfactor == -1:
            turnmsg = "Turn load factor or turn angle not specified in the designbrief dictionary."
            raise ValueError(turnmsg)

        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)

        # W/S at the start of the specified turn test may be less than MTOW/S
        wingloading_pa = wingloading_pa * self.turn_weight_fraction

        twratio, clrequired = self.thrusttoweight_turn(wingloading_pa)

        # What SL T/W will yield the required T/W at the actual altitude?
        temp_c = self.designatm.airtemp_c(self.turnalt_m)
        pressure_pa = self.designatm.airpress_pa(self.turnalt_m)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.turnalt_m)
        turnspeed_mps = self.turnspeed_tas
        mach = self.designatm.mach(turnspeed_mps, self.turnalt_m)
        corr = self._altcorr(temp_c, pressure_pa, mach, density_kgpm3)

        twratio = twratio / corr 

        # Map back to T/MTOW if turn start weight is less than MTOW
        twratio = twratio * self.turn_weight_fraction

        # Which of these points is actually reachable given the clean CLmax?
        feasibletw = np.copy(twratio)
        if self.clmaxclean == -1:
            if self.highlift_clean != -1:
                self.clmaxclean = self.clmax_highlift(self.highlift_clean)
            else:
                clmaxmsg = "CLmaxclean must be specified in the performance dictionary or highlift_type 'clean' in the design dictionary."
                raise ValueError(clmaxmsg)                   
        for idx, val in enumerate(clrequired):
            if val > self.clmaxclean:
                feasibletw[idx] = np.nan

        if len(twratio) == 1:
            return twratio[0], clrequired[0], feasibletw[0]

        return twratio, clrequired, feasibletw


    def twrequired_clm(self, wingloading_pa):
        """Calculates the T/W required for climbing for a range of wing loadings.

        **Parameters**

        wingloading_pa
            float or numpy array, list of wing loading values in Pa.

        **Returns**

        twratio
            array, thrust to weight ratio required for the given wing loadings.

        **See also** ``twrequired``

        **Notes**

        1. Use `twrequired` if a full constraint analysis is desired, as this integrates
        the take-off, turn, climb, cruise, and service ceiling constraints, as well as
        computing the combined constraint boundary.

        2. The calculation currently approximates climb performance on the constant TAS
        assumption (though note that the design brief dictionary variable must specify the
        climb speed as IAS, which is the operationally relevant figure) - a future version
        of the code will remove this approximation and assume constant IAS.

        """

        if self.climbspeed_tas == -1:
            turnmsg = "Climb speed not specified in the designbrief dictionary."
            raise ValueError(turnmsg)

        # Assuming that the climb rate is 'indicated'
        if self.climbrate_mps == -1:
            turnmsg = "Climb rate not specified in the designbrief dictionary."
            raise ValueError(turnmsg)
        climbrate_mps = self.climbrate_mps
        climbrate_mpstroc = climbrate_mps #ORIGINAL_LINE: climbrate_mpstroc = self.designatm.eas2tas(climbrate_mps, self.climbalt_m) !!!!
        
        if self.bpr == -1:
            climbspeed_mpstas = self.bestclimbspeedprop(wingloading_pa, self.climbalt_m)
        else:
            climbspeed_mpstas = self.climbspeed_tas
        
        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)

        # W/S at the start of the specified climb segment may be less than MTOW/S
        wingloading_pa = wingloading_pa * self.climb_weight_fraction

        inddragfact = self.induceddragfact()
        qclimb_pa = self.designatm.dynamicpressure_pa(climbspeed_mpstas, self.climbalt_m)

        cos_sq_theta = (1 - (climbrate_mpstroc / climbspeed_mpstas) ** 2)

        # To be implemented, as 1 + (V/g)*(dV/dh)
        accel_fact = 1.0

        twratio = accel_fact * climbrate_mpstroc / climbspeed_mpstas + \
        (1 / wingloading_pa) * qclimb_pa * self.cdminclean + \
        (inddragfact / qclimb_pa) * wingloading_pa * cos_sq_theta

        # What SL T/W will yield the required T/W at the actual altitude?
        temp_c = self.designatm.airtemp_c(self.climbalt_m)
        pressure_pa = self.designatm.airpress_pa(self.climbalt_m)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.climbalt_m)
        mach = self.designatm.mach(climbspeed_mpstas, self.climbalt_m)
        corr = self._altcorr(temp_c, pressure_pa, mach, density_kgpm3)

        twratio = twratio / corr 

        # Map back to T/MTOW if climb start weight is less than MTOW
        twratio = twratio * self.climb_weight_fraction

        if len(twratio) == 1:
            return twratio[0]

        return twratio , climbspeed_mpstas #ORIGINAL_LINE: return twratio


    def twrequired_sec(self, wingloading_pa):
        """T/W required for a service ceiling for a range of wing loadings"""

        if self.servceil_m == -1:
            secmsg = "Climb rate not specified in the designbrief dictionary."
            raise ValueError(secmsg)

        if self.secclimbspd_tas == -1:
            secmsg = "Best climb speed not specified in the designbrief dictionary."
            raise ValueError(secmsg)
            
        #Calculation of the best climb speed for propeller aircraft
        if self.bpr == -1: 
            secclimbspeed_mpstas = self.bestclimbspeedprop(wingloading_pa, self.servceil_m)
        else:
            secclimbspeed_mpstas = self.secclimbspd_tas

        
        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)

        # W/S at the start of the service ceiling test point may be less than MTOW/S
        wingloading_pa = wingloading_pa * self.sec_weight_fraction
        
        inddragfact = self.induceddragfact()
        qclimb_pa = self.designatm.dynamicpressure_pa(secclimbspeed_mpstas, self.servceil_m)
        
        # Service ceiling typically defined in terms of climb rate (at best climb speed) of
        # dropping to 100feet/min ~ 0.508m/s
        climbrate_mps = co.fpm2mps(100)

        # The climbrate at ceiling is 100 feet/min TAS #ORIGINAL_LINE: # What true climb rate does 100 feet/minute correspond to?
        climbrate_mpstroc = climbrate_mps #ORIGINAL_LINE: climbrate_mpstroc = self.designatm.eas2tas(climbrate_mps, self.servceil_m)
               
        # calculate compressibility drag
        mach = self.designatm.mach(secclimbspeed_mpstas, self.servceil_m)
        cl_wing = wingloading_pa / qclimb_pa        
        cd_c = self.compressibility_drag_wing(mach, cl_wing, self.tc_ratio, self.sweep_025_rad)
        
        # calculate twratio       
        twratio = climbrate_mpstroc / secclimbspeed_mpstas + \
        (1 / wingloading_pa) * qclimb_pa * (self.cdminclean + cd_c) +\
        (inddragfact / qclimb_pa) * wingloading_pa        

        # What SL T/W will yield the required T/W at the actual altitude?
        temp_c = self.designatm.airtemp_c(self.servceil_m)
        pressure_pa = self.designatm.airpress_pa(self.servceil_m)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.servceil_m)
        corr = self._altcorr(temp_c, pressure_pa, mach, density_kgpm3)

        twratio = twratio / corr 

        # Map back to T/MTOW if service ceiling test start weight is less than MTOW
        twratio = twratio * self.sec_weight_fraction

        if len(twratio) == 1:
            return twratio[0]

        return twratio , secclimbspeed_mpstas #ORIGINAL_LINE: return twratio


    def twrequired_crs(self, wingloading_pa):
        """Calculate the T/W required for cruise for a range of wing loadings"""

        if self.cruisespeed_tas == -1:
            cruisemsg = "Cruise speed not specified in the designbrief dictionary."
            raise ValueError(cruisemsg)
        cruisespeed_mps = self.cruisespeed_tas

        if self.cruisealt_m == -1:
            cruisemsg = "Cruise altitude not specified in the designbrief dictionary."
            raise ValueError(cruisemsg)

        wingloading_pa = actools.recastasnpfloatarray(wingloading_pa)

        # W/S at the start of the cruise may be less than MTOW/S
        wingloading_pa = wingloading_pa * self.cruise_weight_fraction

        inddragfact = self.induceddragfact()
        qcruise_pa = self.designatm.dynamicpressure_pa(cruisespeed_mps, self.cruisealt_m)
        
        # calculate compressibility drag
        mach = self.designatm.mach(cruisespeed_mps, self.cruisealt_m)
        cl_wing = wingloading_pa / qcruise_pa
        cd_c = self.compressibility_drag_wing(mach, cl_wing, self.tc_ratio, self.sweep_025_rad)

        # calculate twratio
        twratio = (1 / wingloading_pa) * qcruise_pa * (self.cdminclean + cd_c) + \
        (inddragfact / qcruise_pa) * wingloading_pa

        # What SL T/W will yield the required T/W at the actual altitude?
        temp_c = self.designatm.airtemp_c(self.cruisealt_m)
        pressure_pa = self.designatm.airpress_pa(self.cruisealt_m)
        density_kgpm3 = self.designatm.airdens_kgpm3(self.cruisealt_m)


        corr = self._altcorr(temp_c, pressure_pa, mach, density_kgpm3)

        twratio = twratio / corr 

        # Map back to T/MTOW if cruise start weight is less than MTOW
        twratio = twratio * self.cruise_weight_fraction

        twratio = twratio * (1 / self.cruisethrustfact)

        if len(twratio) == 1:
            return twratio[0]

        return twratio


    def twrequired(self, wingloadinglist_pa, feasibleonly=True):
        """Calculate the T/W required for t/o, trn, clm, crs, sec.

        This method integrates the full set of constraints and it gives the user a
        compact way of performing a full constraint analysis. If a specific constraint
        is required only, the individual methods can be called separately:
        :code:`twrequired_to` (take-off), :code:`twrequired_trn` (turn),
        :code:`twrequired_clm` (climb), :code:`twrequired_trn` (turn),
        :code:`twrequired_crs` (cruise), :code:`twrequired_sec` (service ceiling).

        **Parameters**

        wingloading_pa
            float or numpy array, list of wing loading values in Pa.

        **Returns**

        twreq
            dictionary variable, wherein each entry contains vectors
            related to one of the constraints: :code:`twreq['take-off']`
            (T/W required for take-off), :code:`twreq['liftoffspeed_mps']`
            (liftoff speed in m/s), :code:`twreq['avspeed_mps']` (average
            speed of the take-off run, in m/s), :code:`twreq['turn']`
            (T/W required for the turn), :code:`twreq['turnfeasible']` (same as
            :code:`twreq['turn']`, but with *NaN* where the maximum lift
            coefficient is exceeded), :code:`twreq['turncl']` (lift
            coefficient required in the turn), :code:`twreq['climb']`
            (T/W required for climb), :code:`twreq['cruise']` (T/W required
            for cruise), :code:`twreq['servceil']` (T/W required for the
            service ceiling constraint), :code:`twreq['combined']` (the
            T/W required to meet all of the above).

        """

        tw_to, liftoffspeed_mpstas, avspeed_mpstas = self.twrequired_to(wingloadinglist_pa)
        tw_trn, clrequired, feasibletw_trn = self.twrequired_trn(wingloadinglist_pa)
        tw_clm, climbspeed_mpstas = self.twrequired_clm(wingloadinglist_pa)
        tw_crs = self.twrequired_crs(wingloadinglist_pa)
        tw_sec, secclimbspeed_mpstas = self.twrequired_sec(wingloadinglist_pa)
        if feasibleonly:
            tw_combined = np.amax([tw_to, feasibletw_trn, tw_clm, tw_crs, tw_sec], 0)
        else:
            tw_combined = np.max([tw_to, tw_trn, tw_clm, tw_crs, tw_sec], 0)

        twreq = {
            'take-off': tw_to,
            'liftoffspeed_mpstas': liftoffspeed_mpstas,
            'avspeed_mpstas': avspeed_mpstas,
            'turn': tw_trn,
            'turnfeasible': feasibletw_trn,
            'turncl': clrequired,
            'climb': tw_clm,
            'climbspeed': climbspeed_mpstas,
            'cruise': tw_crs,
            'servceil': tw_sec,
            'secclimbspeed': secclimbspeed_mpstas,
            'combined': tw_combined}

        return twreq


    def powertoweightrequired(self, wingloadinglist_pa, feasibleonly=True):
        """Calculate the power to weight ratio (in HP) required for t/o, trn, clm, crs, sec."""

        twreq = self.twrequired(wingloadinglist_pa, feasibleonly)   #DELETE_AFTER
        if self.bpr != -1 and self.bpr != -2:#DELETE_AFTER
            self.etaprop_to = 1#DELETE_AFTER
            self.etaprop_turn = 1#DELETE_AFTER
            self.etaprop_climb = 1#DELETE_AFTER
            self.etaprop_cruise = 1#DELETE_AFTER
            
        # Take-off power required
        pw_to_wpn = tw2pw(twreq['take-off'], twreq['avspeed_mpstas'], self.etaprop_to)

        # Turn power required
        trnspeed_mpstas = self.turnspeed_tas
        if feasibleonly:
            pw_trn_wpn = tw2pw(twreq['turnfeasible'], trnspeed_mpstas, self.etaprop_turn)
            if np.all(np.isnan(pw_trn_wpn)):
                nanmsg = "All turns are infeasible for the given load factor, speed, and wing loadings."
                warnings.warn(nanmsg, RuntimeWarning)
        else:
            pw_trn_wpn = tw2pw(twreq['turn'], trnspeed_mpstas, self.etaprop_turn)

        # Climb power
        # Conversion to TAS, IAS and EAS conflated, safe for typical prop speeds
        clmspeed_mpstas = twreq['climbspeed']       
        pw_clm_wpn = tw2pw(twreq['climb'], clmspeed_mpstas, self.etaprop_climb)

        # Power for cruise
        crsspeed_mpstas = self.cruisespeed_tas
        pw_crs_wpn = tw2pw(twreq['cruise'], crsspeed_mpstas, self.etaprop_cruise)

        # Power for service ceiling
        # Conversion to TAS, IAS and EAS conflated, safe for typical prop speeds
        secclmspeed_mpstas = twreq['secclimbspeed']
        pw_sec_wpn = tw2pw(twreq['servceil'], secclmspeed_mpstas, self.etaprop_sec)
        
        # Combined power
        pw_combined_wpn = np.amax([pw_to_wpn, pw_trn_wpn, pw_clm_wpn, pw_crs_wpn, pw_sec_wpn], 0)

        pwreq_wpn = {
            'take-off': pw_to_wpn,
            'liftoffspeed_mpstas': twreq['liftoffspeed_mpstas'],
            'avspeed_mpstas': twreq['avspeed_mpstas'],
            'turn': pw_trn_wpn,
            'turncl': twreq['turncl'],
            'climb': pw_clm_wpn,
            'cruise': pw_crs_wpn,
            'servceil': pw_sec_wpn,
            'combined': pw_combined_wpn}

        return pwreq_wpn

   
    def clmax_highlift(self, highlift_type):
        """Calculate the maximum lift coefficient based on the type of high-lift devices and the sweep"""
        
        if highlift_type == None: # No Flaps
            return -0.0002602 * self.sweep_025_deg**2 -0.0008614 * self.sweep_025_deg + 1.51    
            
        elif highlift_type == 'Plain' or highlift_type == 'plain':  # Plain 
            return -0.0002823 * self.sweep_025_deg**2 -0.000141 * self.sweep_025_deg + 1.81  
            
        elif highlift_type == 'Single-slotted' or highlift_type == 'single-slotted': # Single-Slotted 
            return -0.0002599 * self.sweep_025_deg**2 -0.002727 * self.sweep_025_deg + 2.205      
      
        elif highlift_type == 'Fowler' or highlift_type == 'fowler':   # Fowler
            return -0.000283 * self.sweep_025_deg**2 -0.003897 * self.sweep_025_deg + 2.501  
      
        elif highlift_type == 'Double-slotted' or highlift_type == 'double-slotted':    # Double-Slotted            
            return -0.0002574 * self.sweep_025_deg**2 -0.007129 * self.sweep_025_deg + 2.735  
            
        elif highlift_type == 'Double-slotted Fowler' or highlift_type == 'double-slotted Fowler':    # Double-slotted Fowler with Slats
            return -0.0002953 * self.sweep_025_deg**2 -0.006719 * self.sweep_025_deg + 3.014 
            
        elif highlift_type == 'Triple-slotted Fowler' or highlift_type == 'triple-slotted Fowler':    # Triple-slotted Fowler with Slats
            return -0.0003137 * self.sweep_025_deg**2 -0.008903 * self.sweep_025_deg + 3.416 
              
        else:
            highliftmsg = "High-lift device type must be specified in the design dictionary."
            raise ValueError(highliftmsg)
    
    
    def wsmaxlanding_pa(self):
        """Maximum wing loading defined by landing field distance. The method is taken 
        from Scholz (HAW handbook) // Loftin """
        # (W/S)_max = q_vapp * CLmaxApp
        
        if self.clmaxapproach == -1:
            if self.highlift_landing != -1:
                self.clmaxapproach = self.clmax_highlift(self.highlift_landing)
            else:
                clmaxmsg = "CLmaxApp must be specified in the performance dictionary or highlift_type 'landing' in the design dictionary."
                raise ValueError(clmaxmsg)        
        
        if self.landingroll_m == -1:
            landingrollmsg = "Landing roll distance must be specified in the design brief dictionary."
            raise ValueError(landingrollmsg)              
        
        if self.bpr == -1 or self.bpr == -2:
            safety_factor = 1/0.7
        else:
            safety_factor = 1/0.6
            
        landing_field_m = self.landingroll_m * safety_factor

        appspeed_mps = 1.7 * np.sqrt(landing_field_m)
  
        q_pa = self.designatm.dynamicpressure_pa(appspeed_mps, 0)      
        
        return q_pa * self.clmaxapproach
    
    
    def map2static(self):
        """Maps the average take-off thrust to static thrust. If a bypass ratio
        is not specified, it returns a value of 1.
        """
        if self.bpr > 1:
            return (4 / 3) * (4 + self.bpr) / (5 + self.bpr)

        return 1.0
        
    def compressibility_drag_wing(self, mach, cl_wing, tc_ratio, sweep_025_rad): #state,settings,geometry):
        """Computes compressibility drag for a wing

        Assumptions:
        Subsonic to low transonic
        Supercritical airfoil

        Source:
        adg.stanford.edu (Stanford AA241 A/B Course Notes)

        Inputs:
        mach                                                [Unitless]
        cl_wing                                             [Unitless]
        tc_ratio                                            [Unitless]
        sweep_025_rad                                       [radians]

        Outputs:
        cd_c                                                [Unitless]

        Properties Used:
        N/A
        """ 
        
        cos_sweep = np.cos(sweep_025_rad)
        # get effective Cl and sweep
        tc = tc_ratio /(cos_sweep)
        cl = cl_wing / (cos_sweep*cos_sweep)

        # compressibility drag based on regressed fits from AA241
        mcc_cos_ws = 0.922321524499352       \
                   - 1.153885166170620*tc    \
                   - 0.304541067183461*cl    \
                   + 0.332881324404729*tc*tc \
                   + 0.467317361111105*tc*cl \
                   + 0.087490431201549*cl*cl
            
        # crest-critical mach number, corrected for wing sweep
        mcc = mcc_cos_ws / cos_sweep
        # divergence mach number
        MDiv = mcc * ( 1.02 + 0.08*(1 - cos_sweep) )
        
        # divergence ratio
        mo_mc = mach/mcc
        
        # compressibility correlation, Shevell
        dcdc_cos3g = 0.0019*mo_mc**14.641
        
        # compressibility drag
        cd_c = dcdc_cos3g * cos_sweep*cos_sweep*cos_sweep
      
        return cd_c

def tw2pw(thrusttoweight, speed, etap):
    """Converts thrust to weight to power to weight (propeller-driven aircraft)

    **Parameters**

    thrusttoweight
        thrust to weight ratio (non-dimensional)

    speed
        ground speed (in m/s if output in Watts / Newton is required)

    etap
        propeller efficiency (non-dimensional), float

    **Returns**

        power to weight ratio (in W/N if speed is in m/s)

    **See also** ``powerrequired``

    **Notes**

    1. A note on units. If the input speed is in m/s, the other two inputs being
    non-dimensional, the output product is also in m/s, which is equal to W/N
    (W / N = (J/s) / N = (Nm/s) / N = m/s).

    2. The speed input is a kinematic quantity, not an airspeed, so it is generally
    a ground speed (GS) or a true airspeed (TAS) if we are assuming zero wind.

    3. The inputs to the function are scalars or a mix of scalars and `numpy`
    arrays.
    """
    return thrusttoweight * speed / etap

