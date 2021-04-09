## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
# compute_flap_lift.py
#
# Created:  Dec 2013, A. Varyar
# Modified: Feb 2014, T. Orra
#           Jan 2016, E. Botero         

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Units
import numpy as np

# ----------------------------------------------------------------------
#  compute_flap_lift
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_flap_lift(wing):
    """Computes the increase in zero-AoA and max lift due to trailing edge flap deployment

    Assumptions:
    1. The method was proven to work well for wings with moderate AR (between 4 and 16)
    2. Airfoil Cla = 2pi
    3. taper ratio effect on change is CL0 is small, so an area-averaged taper is used
    4. dfClmax = dfCl0
    5. DATCOM CLa uses a quarter-chord sweep as a conservative assumption and assumes M = 0
    6. For CL0, an average value of Cl0 along the span is used

    Source:
    Unknown

    Inputs:
    t_c                 (wing thickness ratio)                 [Unitless]
    flap_type                                                  [string]
    flap_c_chord        (flap chord as fraction of wing chord) [Unitless]
    flap_angle          (flap deflection)                      [radians]
    sweep               (Wing sweep angle)                     [radians]
    wing_Sref           (Wing reference area)                  [m^2]
    wing_affected_area  (Wing area affected by flaps)          [m^2]

    Outputs:
    dcl_max_flaps       (Lift coefficient increase)            [Unitless]

    Properties Used:
    N/A
    """          

    # Unpack and prepare geometry
    AR          = wing.aspect_ratio
    ref_area    = wing.areas.reference
    semispan    = wing.spans.projected*0.5 
    root_chord  = wing.chords.root 

    cf        = wing.flaps.chord     
    bfi       = wing.flaps.span_start 
    bfo       = wing.flaps.span_end  
    flap_type = wing.flaps.type 
    df        = wing.flaps.angle  / Units.degrees

    bf_b = bfo - bfi

    if bf_b != 0:
        if len(wing.Segments.keys())>0: 
            # obtain the geometry for each segment in a loop                                            
            n_segments  = len(wing.Segments.keys())

            chord      = np.zeros(n_segments)
            flap_chord = np.zeros(n_segments)
            loc_span   = np.zeros(n_segments)
            loc_sweep  = np.zeros(n_segments)
            flap_defl  = np.zeros(n_segments)
            Cla_sec    = np.zeros(n_segments)
            Cl0_sec    = np.zeros(n_segments)

            for i in range(n_segments):
                chord[i]     = root_chord * wing.Segments[i].root_chord_percent
                loc_span[i]  = semispan * wing.Segments[i].percent_span_location 
                loc_sweep[i] = wing.Segments[i].sweeps.quarter_chord 

                # Unpack airfoil data
                if wing.airfoils:
                    Cla_sec[i] = wing.Segments[i].Airfoil.airfoil.Cla 
                    Cl0_sec[i] = wing.Segments[i].Airfoil.airfoil.Cl0 
                else: 
                    Cl0_sec[i] = 0.0
                    Cla_sec[i] = 2*np.pi
    
            Cla = np.average(Cla_sec)
            Cl0 = np.average(Cl0_sec)

            # insert flaps into the wing
            cfi       = np.interp(bfi * semispan,loc_span,chord)
            cfo       = np.interp(bfo * semispan,loc_span,chord)
            inb_ind   = loc_span.searchsorted(bfi * semispan)
            loc_span  = np.insert(loc_span, inb_ind, bfi * semispan)
            loc_span  = np.insert(loc_span, inb_ind, bfi * semispan)
            outb_ind  = loc_span.searchsorted(bfo * semispan)
            loc_span  = np.insert(loc_span, outb_ind, bfo * semispan)
            loc_span  = np.insert(loc_span, outb_ind, bfo * semispan)
            chord     = np.insert(chord, inb_ind, cfi)
            chord     = np.insert(chord, inb_ind, cfi)
            chord     = np.insert(chord, outb_ind, cfo)
            chord     = np.insert(chord, outb_ind, cfo)
            flap_chord = np.insert(flap_chord, inb_ind, cf)
            flap_chord = np.insert(flap_chord, inb_ind, 0.0)
            flap_chord = np.insert(flap_chord, outb_ind, cf)
            flap_chord = np.insert(flap_chord, outb_ind+1, 0.0)

            n_segments += 4
            S_segm     = np.zeros(n_segments-1)
            taper      = np.zeros(n_segments-1)
            for i in range(n_segments-1):
                taper[i]  = chord[i+1]/chord[0]
                S_segm[i] = 0.5 * (chord[i]+chord[i+1]) * (loc_span[i+1] - loc_span[i])

        else:
            try:
                chord = np.array([root_chord, root_chord * wing.taper])
            except:
                chord = np.array([root_chord, wing.chords.tip]) 

            # Unpack airfoil data
            Cla = wing.Airfoil.Cla  
            Cl0 = wing.Airfoil.Cl0 

            loc_span  = np.array([0, semispan])
            loc_sweep = np.array([wing.sweeps.quarter_chord, wing.sweeps.quarter_chord])
            taper     = np.array([chord[i+1]/chord[i]])

            cfi       = np.interp(bfi * semispan,loc_span,chord)
            cfo       = np.interp(bfo * semispan,loc_span,chord)
            inb_ind   = loc_span.searchsorted(bfi * semispan)
            loc_span  = np.insert(loc_span, inb_ind, bfi * semispan)
            outb_ind  = loc_span.searchsorted(bfo * semispan)
            loc_span  = np.insert(loc_span, outb_ind, bfo * semispan)
            chord     = np.insert(chord, inb_ind, cfi)
            chord     = np.insert(chord, outb_ind, cfo)

            n_segments = 4
            taper      = np.zeros(n_segments-1)
            for i in range(n_segments-1):
                taper[i]  = chord[i+1]/chord[0]
                S_segm[i] = 0.5 * (chord[i]+chord[i+1]) * (loc_span[i+1] - loc_span[i])

        # insert flaps into the wing
        loc_sweep = np.insert(loc_sweep, inb_ind, loc_sweep[inb_ind-1])
        loc_sweep = np.insert(loc_sweep, inb_ind, loc_sweep[inb_ind-1])
        loc_sweep = np.insert(loc_sweep, outb_ind, loc_sweep[outb_ind-1])
        loc_sweep = np.insert(loc_sweep, outb_ind, loc_sweep[outb_ind-1])
        flap_defl = np.insert(flap_defl, inb_ind, df)
        flap_defl = np.insert(flap_defl, inb_ind, 0)
        flap_defl = np.insert(flap_defl, outb_ind, df)
        flap_defl = np.insert(flap_defl, outb_ind+1, 0)
        flap_defl[inb_ind+1:outb_ind] = df  

        # estimate flap chord extension
        dc_cf = np.zeros(len(flap_defl))
        for i in range(len(flap_defl)):
            if flap_defl[i] != 0:
                if flap_type == 'single_slotted':
                    dc_cf[i] = 5.6096e-3*flap_defl[i]-3.6859e-3
                elif flap_type == 'single_slotted_Fowler' or flap_type == 'double_slotted_fixed_vane':
                    dc_cf[i] = 4.983e-3*flap_defl[i]+4.1435e-1
                elif flap_type == 'double_slotted':
                    dc_cf[i] = 3.5038e-6*np.power(flap_defl[i],3)-4.0912e-4*np.power(flap_defl[i],2)+2.4173e-2*flap_defl[i]-3.3319e-4   
                elif flap_type == 'double_slotted_Fowler' or flap_type == 'triple_slotted_Fowler':
                    dc_cf[i] = 9.111e-3*flap_defl[i]+3.9645e-1
                else:
                    dc_cf[i] = 0.0

        # estimate flap efficiency coefficient
        if flap_type == 'plain' and cf != 0:
            eta = 0.6083 + 1.625*cf + 0.0142*df -2.083*cf**2 -0.0369*cf*df -0.00169*df**2 + 0.05103*cf**2*df - \
                       9.321e-05*cf*df**2 + 4.285e-05*df**3 -0.0002486*cf**2*df**2 + 4.474e-06*cf*df**3 -3.308e-07*df**4    
        elif flap_type == 'single_slotted':
            eta = 3.3938e-8*df**4 - 2.3366e-6*df**3 - 1.2538e-4*df**2 + 7.21e-4*df + 8.1952e-1
        elif flap_type == 'single_slotted_Fowler':
            eta = -8.5591e-6*df**3 + 2.2983e-4*df**2 - 3.1891e-3*df + 8.8817e-1
        elif flap_type == 'double_slotted_fixed_vane':
            eta = 4.9316e-8*df**4 - 6.2302e-6*df**3 + 1.4598e-4*df**2 - 1.2555e-3*df + 7.5178e-1
        elif flap_type == 'triple_slotted':
            eta = -1.4561e-8*df**4 + 1.2379e-6*df**3 - 7.9402e-5*df**2 + 1.1505e-3*df + 8.0942e-1
        else:
            eta = 0

        # Calculate modified chord lengths due to flaps
        chordp = np.multiply(chord,1 + dc_cf * cf)

        # Estimate theoretical flap lift factor
        c_pr_c  = np.divide(chordp,chord)
        theta   = np.arccos(2*np.divide(np.multiply(flap_chord,chord),chordp)-1)
        alphad  = 1 - (theta - np.sin(theta))/np.pi 
        for i in range(len(alphad)):
            if alphad[i] != 0:
                alphad0 = alphad[i]
                c_pr_c0 = c_pr_c[i]
                break

        # Estimate airfoil Cl0 increment due to flaps
        flap_defl = flap_defl * Units.degrees
        dfCl0 = np.multiply(Cla * np.multiply(eta*alphad,flap_defl),c_pr_c) + Cl0 * (c_pr_c-1)

        # average dfCl0, dfClmax, and useful S_segm calculations
        dfCl0_av   = np.zeros(len(dfCl0)-3)
        S_segm_av  = np.zeros(len(dfCl0)-3)
        sweep_segm = np.zeros(len(dfCl0)-3)
        taper_segm = np.zeros(len(dfCl0)-3)
        ind = 0
        for i in range(len(dfCl0_av)):
            dfCl0_av[i]   = np.average([dfCl0[ind],dfCl0[ind+1]])*np.average(np.cos([loc_sweep[ind],loc_sweep[ind+1]]))
            sweep_segm[i] = np.average([loc_sweep[ind],loc_sweep[ind+1]])
            ind +=2
        ind = 0
        for i in range(len(S_segm)):
            if S_segm[i] != 0:
                S_segm_av[ind]  = S_segm[i]
                taper_segm[ind] = taper[i]
                ind += 1

        ## Compute increase in CLmax
        dfCLmax = 0.92 * np.sum(np.multiply(dfCl0_av,S_segm_av))/(0.5 * ref_area) 

        ## Compute increase in CL0 
        # calculate Cl0 for any section
        for i in range(len(dfCl0_av)):
            if dfCl0_av[i] != 0:
                dfCl0av = dfCl0_av[i]
                break

        # Estimate area-averaged quarter-chord sweep
        sweep_av  = np.sum(np.multiply(sweep_segm, S_segm_av)) / np.sum(S_segm_av)
        sum_total = np.sum(S_segm_av) 

        # Calculate the wing lift-curve-slope
        kappa = Cla / (2*np.pi)
        CLa   = 2*np.pi*AR / (2 + (4 + (AR/kappa)**2*(1+(np.tan(sweep_av))**2))**0.5)

        # Calculate an area-averaged taper ratio
        taper_av = np.sum(np.multiply(taper_segm, S_segm_av)) / np.sum(S_segm_av)

        # flap span ratio correction
        Kb = 0.0001109 + 1.236*bf_b + 0.06372*taper_av + 0.145*bf_b**2 + 0.02106*bf_b*taper_av \
                   -0.06372*taper_av**2 -0.38*bf_b**3 + 7.75e-16*bf_b**2*taper_av -0.02106*bf_b*taper_av**2

        # AR correction
        aCl = alphad0 * eta
        Kc  = 2.041 -0.1186*AR -3.964*aCl + 0.005525*AR**2 + 0.336*AR*aCl + 7.059*aCl**2 -9.937e-05*AR**3  \
                -0.009523*AR**2*aCl -0.3918*AR*aCl**2 -6.292*aCl**3 + 6.576e-05*AR**3*aCl  \
                +0.005236*AR**2*aCl**2 + 0.1699*AR*aCl**3 + 2.189*aCl**4

        # Final CL increment due to flap
        dfCL0 = dfCl0av * CLa/(2*np.pi) * Kc * Kb     

        # Pack output
        wing.sweeps.quarter_chord = sweep_av

    else:
        dfCL0     = 0.0
        dfCLmax   = 0.0
        alphad0   = 0.0
        c_pr_c0   = 0.0
        sum_total = 0.0
        Cla       = 2 * np.pi

    return dfCLmax, dfCL0, alphad0, c_pr_c0, sum_total, Cla


# ----------------------------------------------------------------------
#   Module Tests
# ----------------------------------------------------------------------
# this will run from command line, put simple tests for your code here
if __name__ == '__main__':

    #imports
    import SUAVE
    from SUAVE.Core import Units

    # Test case
    t_c             = 0.11
    flap_type       = 'single_slotted'
    flap_chord      = 0.28
    flap_angle      = 30. * Units.deg
    sweep           = 30. * Units.deg
    wing_Sref       = 120.
    wing_flap_area  = 120. * .6

    dcl_flap = compute_flap_lift(t_c,flap_type,flap_chord,flap_angle,sweep,wing_Sref,wing_flap_area)
    print('Delta CL due to Flaps: ', dcl_flap)

    # Test case
    t_c             = 0.11
    flap_type       = 'none'
    flap_chord      = 0.
    flap_angle      = 0. * Units.deg
    sweep           = 25. * Units.deg
    wing_Sref       = 120.
    wing_flap_area  = 0.

    dcl_flap = compute_flap_lift(t_c,flap_type,flap_chord,flap_angle,sweep,wing_Sref,wing_flap_area)
    print('Delta CL due to Flaps: ', dcl_flap)