## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
# compute_slat_lift.py
#
# Created:  Dec 2013, A. Variyar
# Modified: Feb 2014, T. Orra
#           Jun 2014, T. Orra 
#           Jan 2016, E. Botero

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import SUAVE
from SUAVE.Core import Units
import numpy as np

# ----------------------------------------------------------------------
#  compute_slat_lift
# ----------------------------------------------------------------------

## @ingroup Methods-Aerodynamics-Fidelity_Zero-Lift
def compute_slat_lift(wing,w_Clmax):
    """Computes the increase in zero-AoA and max lift due to leading_edge devices deployment

    Assumptions:

    Source:
    adg.stanford.edu (Stanford AA241 A/B Course Notes)

    Inputs:
    slat_angle   [radians]
    sweep_angle  [radians]

    Outputs:
    dcl_slat     [Unitless]

    Properties Used:
    N/A
    """   

    # Unpack and prepare geometry
    ref_area    = wing.areas.reference
    semispan    = wing.spans.projected*0.5 
    root_chord  = wing.chords.root 

    cs        = wing.slats.chord     
    bsi       = wing.slats.span_start 
    bso       = wing.slats.span_end  
    ds        = wing.slats.angle 

    if bso - bsi != 0:
        if len(wing.Segments.keys())>0: 
            # obtain the geometry for each segment in a loop                                            
            n_segments  = len(wing.Segments.keys())

            chord      = np.zeros(n_segments)
            slat_chord = np.zeros(n_segments)
            loc_span   = np.zeros(n_segments)
            loc_sweep  = np.zeros(n_segments)
            slat_defl  = np.zeros(n_segments)
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
            csi       = np.interp(bsi * semispan,loc_span,chord)
            cso       = np.interp(bso * semispan,loc_span,chord)
            inb_ind   = loc_span.searchsorted(bsi * semispan)
            loc_span  = np.insert(loc_span, inb_ind, bsi * semispan)
            loc_span  = np.insert(loc_span, inb_ind, bsi * semispan)
            outb_ind  = loc_span.searchsorted(bso * semispan)
            loc_span  = np.insert(loc_span, outb_ind, bso * semispan)
            loc_span  = np.insert(loc_span, outb_ind, bso * semispan)
            chord     = np.insert(chord, inb_ind, csi)
            chord     = np.insert(chord, inb_ind, csi)
            chord     = np.insert(chord, outb_ind, cso)
            chord     = np.insert(chord, outb_ind, cso)
            slat_chord = np.insert(slat_chord, inb_ind, cs)
            slat_chord = np.insert(slat_chord, inb_ind, 0.0)
            slat_chord = np.insert(slat_chord, outb_ind, cs)
            slat_chord = np.insert(slat_chord, outb_ind+1, 0.0)

            n_segments += 4

            S_segm    = np.zeros(n_segments-1)
            for i in range(n_segments-1):
                S_segm[i] = 0.5 * (chord[i]+chord[i+1]) * (loc_span[i+1] - loc_span[i])

        else:
            try:
                chord = np.array([root_chord, root_chord * wing.taper])
            except:
                chord = np.array([root_chord, wing.chords.tip]) 

            # Unpack airfoil data
            try:
                Cla = wing.Airfoil.Cla  
                Cl0 = wing.Airfoil.Cl0 
            except:
                Cla = 2 * np.pi
                Cl0 = 0.0

            loc_span  = np.array([0, semispan])
            loc_sweep = np.array([wing.sweeps.quarter_chord, wing.sweeps.quarter_chord])

            # insert flaps into the wing
            slat_chord = np.zeros(2)
            csi        = np.interp(bsi * semispan,loc_span,chord)
            cso        = np.interp(bso * semispan,loc_span,chord)
            inb_ind    = loc_span.searchsorted(bsi * semispan)
            loc_span   = np.insert(loc_span, inb_ind, bsi * semispan)
            loc_span   = np.insert(loc_span, inb_ind, bsi * semispan)
            outb_ind   = loc_span.searchsorted(bso * semispan)
            loc_span   = np.insert(loc_span, outb_ind, bso * semispan)
            loc_span   = np.insert(loc_span, outb_ind, bso * semispan)
            chord      = np.insert(chord, inb_ind, csi)
            chord      = np.insert(chord, inb_ind, csi)
            chord      = np.insert(chord, outb_ind, cso)
            chord      = np.insert(chord, outb_ind, cso)
            slat_chord = np.insert(slat_chord, inb_ind, cs)
            slat_chord = np.insert(slat_chord, inb_ind, 0.0)
            slat_chord = np.insert(slat_chord, outb_ind, cs)
            slat_chord = np.insert(slat_chord, outb_ind+1, 0.0)

            n_segments = 6

            S_segm    = np.zeros(n_segments-1)
            for i in range(n_segments-1):
                S_segm[i] = 0.5 * (chord[i]+chord[i+1]) * (loc_span[i+1] - loc_span[i])

        # insert sweep and slat deflections into the wing
        loc_sweep = np.insert(loc_sweep, inb_ind, loc_sweep[inb_ind-1])
        loc_sweep = np.insert(loc_sweep, inb_ind, loc_sweep[inb_ind-1])
        loc_sweep = np.insert(loc_sweep, outb_ind, loc_sweep[outb_ind-1])
        loc_sweep = np.insert(loc_sweep, outb_ind, loc_sweep[outb_ind-1])
        slat_defl = np.insert(slat_defl, inb_ind, ds)
        slat_defl = np.insert(slat_defl, inb_ind, 0)
        slat_defl = np.insert(slat_defl, outb_ind, ds)
        slat_defl = np.insert(slat_defl, outb_ind+1, 0)
        slat_defl[inb_ind+1:outb_ind] = ds 

        # Estimate theoretical flap lift factor
        theta   = np.arccos(1-2*slat_chord)
        dsCl0   = np.multiply(-(theta - np.sin(theta))/np.pi,slat_chord) * 2 * np.pi
        # average dsCl0 and useful S_segm calculations
        dsCl0_av  = np.zeros(len(dsCl0)-3)
        S_segm_av = np.zeros(len(dsCl0)-3)
        ind = 0
        for i in range(len(dsCl0_av)):
            dsCl0_av[i] = np.average([dsCl0[ind],dsCl0[ind+1]])*np.average(np.cos([loc_sweep[ind],loc_sweep[ind+1]]))
            ind +=2
        ind = 0
        for i in range(len(S_segm)):
            if S_segm[i] != 0:
                S_segm_av[ind] = S_segm[i]
                ind += 1

        dsClmax = np.zeros(len(dsCl0_av))
        for i in range(len(dsCl0_av)):
            if dsCl0_av[i] != 0:
                dsClmax[i] = 0.93 * (Cl0+dsCl0_av[i]+0.47*Cla)/(1+0.035*Cla) - w_Clmax
            
        dsCL0   = 0.92 * np.sum(np.multiply(dsCl0_av,S_segm_av)) / (0.5 * ref_area)
        dsCLmax = 0.92 * np.sum(np.multiply(dsClmax,S_segm_av)) / (0.5 * ref_area) 

        S_slat_tot = np.sum(S_segm_av) 

    else:
        dsCL0      = 0.0
        dsCLmax    = 0.0
        S_slat_tot = 0.0

    return dsCL0, dsCLmax, S_slat_tot      

# ----------------------------------------------------------------------
#   Module Tests
# ----------------------------------------------------------------------
# this will run from command line, put simple tests for your code here
if __name__ == '__main__':

    #imports
    import pylab as plt
    import matplotlib
    matplotlib.interactive(True)
    import scipy as sp
    import SUAVE
    from SUAVE.Core import Units

    #define array of sweep and deflection
    sweep_vec = sp.linspace(-10,30,20) * Units.deg
    deflection_vec = sp.linspace(0,50,6)* Units.deg

    dcl_slat = sp.zeros((len(sweep_vec),len(deflection_vec)))
    legend = ''

    for i in range(len(sweep_vec)):
        for j in range(len(deflection_vec)):
            sweep = sweep_vec[i]
            deflection = deflection_vec[j]
            dcl_slat[i,j] = compute_slat_lift(deflection,sweep)

    # ------------------------------------------------------------------
    #   Plotting Delta CL due to Slat vs Sweep angle
    # ------------------------------------------------------------------
    title = "Delta dCL_slat vs Wing sweep"
    plt.figure(1); 
    for deflection in range(len(deflection_vec)):
        plt.plot(sweep_vec/Units.deg , dcl_slat[:,deflection] ,'bo-', \
                    label = 'Deflection: ' +  str(deflection_vec[deflection]/Units.deg) + ' deg')
    plt.xlabel('Sweep angle (deg)'); plt.ylabel('delta CL due to Slat')
    plt.title(title); plt.grid(True)
    legend = plt.legend(loc='upper right', shadow = 'true')
    plt.show(block=True)

