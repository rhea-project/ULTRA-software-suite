## @ingroup methods-Aerodynamics-Quasi_3D
# Quasi_3D_wing.py
#
# Created: Jul 2019, S. Karpuk 

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Units, Data
from SUAVE.Core import redirect

from SUAVE.Analyses import Process
from SUAVE.Analyses.Mission.Segments.Conditions.Aerodynamics                   import Aerodynamics
from SUAVE.Analyses.Mission.Segments.Conditions.Conditions                     import Conditions
from SUAVE.Methods.Geometry.Three_Dimensional.compute_wing_section_coordinates import compute_wing_section_coordinates as wing_sections
from SUAVE.Methods.Geometry.Three_Dimensional.Compute_sweep                    import Compute_sweep
from SUAVE.Methods.Aerodynamics.AVL.Data.Settings                              import Settings
from SUAVE.Methods.Aerodynamics.AVL.read_lift_distribution                     import read_lift_distribution
from SUAVE.Methods.Aerodynamics.Xfoil.Data.Settings                            import Settings as Airfoil_Settings
from SUAVE.Methods.Aerodynamics.Quasi_3D.Data.Settings                         import Settings as Q3D_Settings
from SUAVE.Methods.Aerodynamics.Xfoil.initialize_inputs                        import initialize_inputs
from SUAVE.Methods.Aerodynamics.Xfoil.write_input_deck                         import write_input_deck
from SUAVE.Methods.Aerodynamics.Xfoil.run_analysis                             import run_analysis
from SUAVE.Methods.Aerodynamics.Xfoil.read_results                             import read_results

# Package imports
import os
import glob
import math
import time
import shutil
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

def Quasi_3D_wing(self):
    """Computes the parasite drag due to wings using a Quasi-3D approach 

    Assumptions:
    Basic fit

    Source:
    Mariens, Elham, Tooren, "Quasi-Three-Dimensional Aerodynamic Solver for
    Multidisciplinary Design Optimization of Lifting Surfaces",
    Journal of Aircraft, Vol 51, #2, March-April 2014

    Inputs:
    settings.wing_parasite_drag_form_factor      [Unitless]
    state.conditions.freestream.
      mach_number                                [Unitless]
      temperature                                [K]
      reynolds_number                            [Unitless]
    geometry.
      areas.reference                            [m^2]
      chords.mean_aerodynamic                    [m]
      thickness_to_chord                         [Unitless]
      sweeps.quarter_chord                       [radians]
      aspect_ratio                               [Unitless]
      spans.projected                            [m]
      areas.exposed                              [m^2]
      areas.affected                             [m^2]
      areas.wetted                               [m^2]
      transition_x_upper                         [Unitless]
      transition_x_lower                         [Unitless]
      
      
    Outputs:
    CD                                           [Unitless]

    Properties Used:
    N/A
    
    """
    
    # Unpack
    wing           = self.geometry.wings.main_wing
    settings       = Settings()
    drag_sect      = self.number_of_drag_setions
    file_tag       = self.geometry._base.tag
    ref_line       = self.ref_line_percent
      
    airfoil_settings   = Airfoil_Settings()
    Q3D_settings       = Q3D_Settings()
    run_folder         = os.path.abspath(settings.filenames.run_folder)
    run_Q3Dfolder      = os.path.abspath(Q3D_settings.filenames.run_folder)
    run_airfoil_folder = os.path.abspath(airfoil_settings.filenames.run_folder)
    geometry_file      = file_tag + '.avl'

    AoA  = self.training.angle_of_attack        
    Mach = self.training.Mach                   
    Re   = self.training.Re

    CD         = np.zeros([len(AoA)*len(Mach)*len(Re),1])
    table_size = len(AoA)*len(Mach)*len(Re)

    # Obtain Wing airfoil section numpy array from the geometry file
    #-----------------------------------------------------------------------------------
    with redirect.folder(run_folder,force=False):
        wing_sections(geometry_file,wing)
    
    n_segments  = len(wing.Segments.keys())
    Strip_points, Y_strip, Trans_sweep, chord, tilted_chord, drag_sect = prepare_airfoils(n_segments,wing,file_tag,drag_sect,ref_line)

    # Copy all AVL files from the AVL folder to the Q3D folder
    #-----------------------------------------------------------------------------------
    files = os.listdir(run_folder)
    try:   
        os.mkdir(run_Q3Dfolder)
    except:
        shutil.rmtree(run_Q3Dfolder)
        os.mkdir(run_Q3Dfolder)
        
    for f in files:
        shutil.copy(run_folder + '\\' + f, run_Q3Dfolder)

    
    # Interpolate the lift distribution
    #-----------------------------------------------------------------------------------
    with redirect.folder(run_Q3Dfolder,force=False):
        Lift_distribution = read_lift_distribution(AoA,Mach,Y_strip)
        
    # Delete all AVL run files
    #-----------------------------------------------------------------------------------
    filenames = ['commands', 'batch', 'results', 'Trefftz']
    with redirect.folder(run_Q3Dfolder,force=False):
        for iii in range(len(filenames)):  
            for rem_filename in glob.glob(filenames[iii] + '*'):
                print(rem_filename)
                os.remove(rem_filename)
        
        
    # Calculate local sweep for each strip
    #-----------------------------------------------------------------------------------
    Lift_distribution = compute_strip_sweep(wing,Trans_sweep,Lift_distribution)

    # Evaluate section drag using Q3D
    #-----------------------------------------------------------------------------------
    time0 = time.time()
    print('Running Q3D Analysis: \nTotal number of cases = ' + str(table_size * drag_sect) + '...')

    time0 = time.time()
    for ii in range(len(Lift_distribution.keys())):
        CD = Calculate_Drag(self,ii,Lift_distribution,tilted_chord,chord,Strip_points,ref_line,wing, \
                            CD,Re,drag_sect,run_airfoil_folder,Y_strip)
                 
    time1 = time.time()
    print('Total time = ' + str(time1-time0))
    

    return CD

def Calculate_Drag(self,ii,Lift_distribution,tilted_chord,chord,Strip_points,ref_line,wing,CD,Re,drag_sect,run_airfoil_folder,Y_strip):
    """Runs the Q3D sweep for each free-stream condition

    Assumptions:

    Source:

    Inputs:
    Lift_distribution                           [Data structure] includes information regarrding lift distribution
    tilted_chord                                [m]
    chord                                       [m]
    Strip_points                                [m]
    ref_line                                    [Unitless]
    wing
    CD                                          [Unitless]
    Re                                          [Unitless]
    drag_sect                                   [Unitless]
    run_airfoil_folder                          String
    Y_strip                                     [m]
      
    Outputs:
    CD                                          [Unitless]

    Properties Used:
    N/A
    
    """
    EPS   = 10E-4
    for jj in range(len(Re)):
        Cd_c = np.zeros(drag_sect)
        print('Mach= ' + str(Lift_distribution[ii].Mach) + ' AoA= ' + str(Lift_distribution[ii].AoA) \
                  + ' Re= ' + str(Re[jj]))
        for kk in range(drag_sect):
            alpha   = Lift_distribution[ii].AoA * Units.degrees
            sweep   = Lift_distribution[ii].sweep[kk]
            Mp      = Lift_distribution[ii].Mach * math.cos(sweep)
            Clp     = Lift_distribution[ii].ClQ3D[kk] / (math.cos(sweep))**2
            
            alpha_i       = 0
            Cd_eff        = 0
            error         = 1
            iteration     = 1
            iteration_max = 25 
            while abs(error) > EPS and iteration <= iteration_max:
                M_eff    = Mp / math.cos(alpha_i)
                Cl_eff   = Clp * (math.cos(alpha_i))**2 + Cd_eff * math.sin(alpha_i) / math.cos(alpha_i)
                Re_eff   = Re[jj] * math.cos(sweep) / math.cos(alpha_i) * tilted_chord[kk+1]
                
                diverg_flag  = False
                # Initialize Xfoil input files
                with redirect.folder(run_airfoil_folder,force=False):
                    xfoil_inputs = initialize_inputs(Strip_points,M_eff,Cl_eff,Re_eff,ref_line,wing)

                    # RUN Xfoil!
                    # try different number of airfoil nodes
                    N_panels    = 180
                    N_panel_max = 260
                    er_pan = 1
                    while er_pan > 0 and (N_panels <= N_panel_max):
                        N_panels_old = N_panels
                        xfoil_inputs = write_input_deck(xfoil_inputs,xfoil_inputs.airfoils.Segments[kk],M_eff,Cl_eff,Re_eff,diverg_flag,N_panels)
                        run_analysis(self,xfoil_inputs.airfoils.Segments[kk].result_file, \
                                                                   xfoil_inputs.airfoils.Segments[kk].command_file)
                        results_xfoil, diverg_flag, N_panels = read_results(xfoil_inputs.airfoils.Segments[kk].result_file,diverg_flag,  \
                                                                            Cl_eff,N_panels,N_panel_max)
                        er_pan = N_panels - N_panels_old

                    N_panels = 200
                    if diverg_flag == True:
                        xfoil_inputs = write_input_deck(xfoil_inputs,xfoil_inputs.airfoils.Segments[kk],M_eff,Cl_eff,Re_eff,diverg_flag,N_panels)
                        run_analysis(self,xfoil_inputs.airfoils.Segments[kk].result_file, \
                                                               xfoil_inputs.airfoils.Segments[kk].command_file)
                        results_xfoil, diverg_flag, N_panels = read_results(xfoil_inputs.airfoils.Segments[kk].result_file,diverg_flag,  \
                                                                            Cl_eff,N_panels,N_panel_max)
                    
                    if diverg_flag == True:
                        raise Exception('Cl = ' + str(Cl_eff) + '\n The strip ' + str(kk) + ' airfoil has reached its Clmax. Change the alpha-sweep or the planform shape')
             
                    incidence   = xfoil_inputs.airfoils.Segments[kk].twist
                    alpha_eff   = results_xfoil.alpha_eff
                    Cd_eff      = results_xfoil.Cd_eff
                    Cdpr_eff    = results_xfoil.Cdp_eff
                    alpha_i_old = alpha_i
                alpha_i   = -alpha_eff + (alpha + incidence)/math.cos(sweep)
                error     = alpha_i - alpha_i_old
                iteration = iteration + 1

            Cdf      = (Cd_eff - Cdpr_eff) / (math.cos(alpha_i))        #CHECK THIS FORMULATION
            Cdpr     = Cd_eff * (math.cos(sweep))**3/ math.cos(alpha_i)
            Cd       = Cdf + Cdpr
            Cd_c[kk] = Cd * chord[kk+1]

        CD[ii*len(Re)+jj] = 2/wing.areas.reference * np.trapz(Cd_c,x = Y_strip[1:len(Y_strip)-1])

    return CD

            
def prepare_airfoils(n_segments,wing,file_tag,drag_sect,ref_line):
    """Prepares airfoil sections along the wing perpendicular to the sweep line

    Assumptions:

    Source:

    Inputs:
    wing
    n_segments                                   [Unitless]
    file_tag                                     string
    drag_sect                                    [Unitless]
    ref_line                                     [Unitless]
      
    Outputs:
    Strip_points                                 [m]
    Y_strip                                      [m]
    Trans_sweep                                  [degrees]
    c                                            [m]
    tilted_chord                                 [m]

    Properties Used:
    N/A
    
    """
    b                  = np.zeros(n_segments)
    Trans_sweep        = np.zeros(n_segments)
    LE_sweep           = np.zeros(n_segments)
    Ct                 = np.zeros(n_segments)
    Cr                 = wing.chords.root
    Ct[0]              = Cr
    dihedral           = []
    sweep_for_dihidral = []
    dihedral_span      = [0.0]
    n_strips           = np.zeros(n_segments)
    dY_strip           = np.zeros(drag_sect)

    half_span     = 0.5*wing.spans.projected
    # Calculate local sweeps and span lengths for each segment
    #--------------------------------------------------------------------------------------------------------------------------------------------
    for i_segs in range(1,n_segments):
        Ct[i_segs] = Cr * wing.Segments[i_segs].root_chord_percent
        # Specify local segment span
        b[i_segs] = 0.5 * (wing.Segments[i_segs].percent_span_location - wing.Segments[i_segs-1].percent_span_location) * wing.spans.projected
             
        # Calculate sweeps for each segment
        if wing.Segments[i_segs].sweeps.leading_edge > 0:
            LE_sweep[i_segs]     = wing.Segments[i_segs-1].sweeps.leading_edge
            Trans_sweep[i_segs]  = Compute_sweep(b[i_segs],Ct[i_segs-1],Ct[i_segs],ref_line,LE_sweep[i_segs] ,'from leading')
        else:
            Quart_sweep = wing.Segments[i_segs-1].sweeps.quarter_chord
            LE_sweep[i_segs]     = Compute_sweep(b[i_segs],Ct[i_segs-1],Ct[i_segs],0.25,Quart_sweep,'to leading')      # Find the LE sweep from the quarter-chord sweep
            if ref_line == 0.5:
                Trans_sweep[i_segs]  = Compute_sweep(b[i_segs],Ct[i_segs-1],Ct[i_segs],ref_line,LE_sweep[i_segs] ,'from leading')
            else:
                Trans_sweep[i_segs]  = Quart_sweep

    # Discretize the wing into strips
    #--------------------------------------------------------------------------------------------------------------------------------------------
    if wing.manual_segments == True:
        drag_sect = 0
        for i in range(n_segments-1):
            drag_sect = drag_sect + wing.Segments[i].sections_outboard
        Y_strip = []

        for i in range(n_segments-1):
            b1       = 0.5 * wing.Segments[i].percent_span_location * wing.spans.projected
            b2       = 0.5 * wing.Segments[i+1].percent_span_location * wing.spans.projected
            dy       = (b2-b1) / (wing.Segments[i].sections_outboard+1)
            Y_strip1 = np.arange(b1,b2,dy)

            # Check min spacing violation for adjacent strips
            if Trans_sweep[i] > 0:
                h0 = 0.375*Ct[i]*abs(math.sin(2*Trans_sweep[i]))/math.cos(wing.Segments[i].dihedral_outboard) \
                     + 0.5 * wing.spans.projected * 0.01
            else:
                h0 = 0.125*Ct[i]*abs(math.sin(2*Trans_sweep[1]))/math.cos(wing.Segments[0].dihedral_outboard) + 0.5 * wing.spans.projected * 0.01

            if dy < h0:
                raise Exception('Strips in Section ' + str(i) + ' intersect the section boundaries. Change the number of segments there')
            
            if i == 0:
                Y_strip = np.append(Y_strip,Y_strip1)
            elif i == n_segments-2:
                Y_strip = np.append(Y_strip,Y_strip1[1:])
            else:
                Y_strip = np.append(Y_strip,Y_strip1[1:])
        Y_strip = np.append(Y_strip,0.5 * wing.spans.projected * wing.Segments[n_segments-1].percent_span_location)

    else:
        Y_strip = np.zeros(drag_sect)
        dy      = 0.5 * wing.spans.projected * wing.Segments[n_segments-1].percent_span_location / (drag_sect + 1)
        Y_strip = np.arange(0,0.5*wing.spans.projected,dy)

        # check the tip and root spacing to fit all strips completely
        if Trans_sweep[1] > 0:
            h0 = 0.375*Cr*abs(math.sin(2*Trans_sweep[1]))/math.cos(wing.Segments[0].dihedral_outboard) \
                 + 0.5 * wing.spans.projected * 0.01
            if Y_strip[1] < h0:
                Y_strip[1] = h0
                dy         = (0.5 * wing.spans.projected * wing.Segments[n_segments-1].percent_span_location - h0)/ (drag_sect)
                for j in range(2,len(Y_strip)-1):
                    Y_strip[j] = Y_strip[j-1] + dy
        else:
            h0 = 0.125*Cr*abs(math.sin(2*Trans_sweep[1]))/math.cos(wing.Segments[0].dihedral_outboard) + 0.5 * wing.spans.projected * 0.01
            if Y_strip[1] < h0:
                Y_strip[1] = h0
                dy         = (0.5 * wing.spans.projected * wing.Segments[n_segments-1].percent_span_location - h0)/ (drag_sect)
                for j in range(2,len(Y_strip)-1):
                    Y_strip[j] = Y_strip[j-1] + dy
            
        if Trans_sweep[n_segments-1] > 0:        
            h0 = 0.125*Ct[n_segments-1]*abs(math.sin(2*Trans_sweep[n_segments-1]))/math.cos(wing.Segments[0].dihedral_outboard) \
                 + 0.5 * wing.spans.projected * 0.01
            if dy < h0:
                Y_strip[len(Y_strip)-1] = 0.5 * wing.spans.projected - h0
        else:
            h0 = 0.375*Ct[n_segments-1]*abs(math.sin(2*Trans_sweep[n_segments-1]))/math.cos(wing.Segments[0].dihedral_outboard)  \
                 + 0.5 * wing.spans.projected * 0.01

            if dy < h0:
                Y_strip[len(Y_strip)-1] = 0.5 * wing.spans.projected - h0
  
    # Create an equation of the quarter/half chord, leading and trailing edge lines
    #--------------------------------------------------------------------------------------------------------------------------------------------
    K_LE     = np.zeros([n_segments,2])
    K_Trans  = np.zeros([n_segments,2])
    K_TE     = np.zeros([n_segments,2])
    N_points = np.shape(wing.Segments[0].Airfoil.airfoil.points)
 
    for i in range(1,n_segments):
        X_trans   = wing.Segments[i].origin[0] + ref_line*(wing.Segments[i].Airfoil.airfoil.points[0][0]-wing.Segments[i].origin[0])
        Y_trans   = wing.Segments[i].origin[1]
        X_trans_1 = wing.Segments[i-1].origin[0] + ref_line*(wing.Segments[i-1].Airfoil.airfoil.points[0][0]-wing.Segments[i-1].origin[0])
        Y_trans_1 = wing.Segments[i-1].origin[1]

        K_LE[i][0]    = (wing.Segments[i-1].origin[0] - wing.Segments[i].origin[0])/  \
                        (wing.Segments[i-1].origin[1] - wing.Segments[i].origin[1])
        K_TE[i][0]    = (wing.Segments[i-1].Airfoil.airfoil.points[0][0] - wing.Segments[i].Airfoil.airfoil.points[0][0])/  \
                        (wing.Segments[i-1].Airfoil.airfoil.points[0][1] - wing.Segments[i].Airfoil.airfoil.points[0][1])
        K_Trans[i][0] = (X_trans_1-X_trans)/(Y_trans_1-Y_trans)
        K_LE[i][1]    = (wing.Segments[i].origin[1]*wing.Segments[i-1].origin[0]-  \
                         wing.Segments[i-1].origin[1]*wing.Segments[i].origin[0])/   \
                        ((-wing.Segments[i-1].origin[1] + wing.Segments[i].origin[1]))
        K_TE[i][1]    = (wing.Segments[i].Airfoil.airfoil.points[0][1]*wing.Segments[i-1].Airfoil.airfoil.points[0][0]-  \
                         wing.Segments[i-1].Airfoil.airfoil.points[0][1]*wing.Segments[i].Airfoil.airfoil.points[0][0])/   \
                        ((-wing.Segments[i-1].Airfoil.airfoil.points[0][1] + wing.Segments[i].Airfoil.airfoil.points[0][1]))
        K_Trans[i][1] = (-Y_trans*X_trans_1+X_trans*Y_trans_1)/(-Y_trans+Y_trans_1)       

    # Create line equations for each chord strip
    #--------------------------------------------------------------------------------------------------------------------------------------------
    K               = np.zeros([drag_sect+1,2])
    tilted_chord_LE = np.zeros([drag_sect+1,2])
    tilted_chord_TE = np.zeros([drag_sect+1,2])
    tilted_chord    = np.zeros(drag_sect+1)
    c               = np.zeros(drag_sect+1)

    for i in range(1,drag_sect+1):
        sect_index = 1
        while Y_strip[i] > wing.Segments[sect_index].origin[1]:
            sect_index = sect_index + 1
        K[i][0] = -1/K_Trans[sect_index][0]
        c[i]    = (K_TE[sect_index][0]*Y_strip[i] + K_TE[sect_index][1]) - (K_LE[sect_index][0]*Y_strip[i] + K_LE[sect_index][1])
        X_trans = (K_LE[sect_index][0]*Y_strip[i] + K_LE[sect_index][1]) + ref_line * c[i]
        K[i][1] = X_trans - Y_strip[i]*K[i][0]

        # Find leading and trailing edges of each strip (with the neighbouring section check)
        #-------------------------------------------------------------------------------------------
        # Check the special case - rectangular wing with zero sweep
        if K_Trans[sect_index][0] == 0 and K_LE[sect_index][0] == 0 and K_TE[sect_index][0] == 0:
            tilted_chord_LE[i][0] = K_LE[sect_index][1]
            tilted_chord_LE[i][1] = Y_strip[i]
            tilted_chord_TE[i][0] = K_TE[sect_index][1]
            tilted_chord_TE[i][1] = Y_strip[i]
        else:
            if Trans_sweep[sect_index] > 0:
                tilted_chord_LE[i][1] = -(K[i][1] - K_LE[sect_index][1]) / (K[i][0] - K_LE[sect_index][0])
                tilted_chord_LE[i][0] = K[i][0] * tilted_chord_LE[i][1] + K[i][1]
                if tilted_chord_LE[i][1] > wing.Segments[sect_index].origin[1]:
                    tilted_chord_LE[i][1] = -(K[i][1] - K_LE[sect_index+1][1]) / (K[i][0] - K_LE[sect_index+1][0])
                    tilted_chord_LE[i][0] = K[i][0] * tilted_chord_LE[i][1] + K[i][1]

                tilted_chord_TE[i][1] = -(K[i][1] - K_TE[sect_index][1]) / (K[i][0] - K_TE[sect_index][0])
                tilted_chord_TE[i][0] = K[i][0] * tilted_chord_TE[i][1] + K[i][1]  
                if tilted_chord_TE[i][1] < wing.Segments[sect_index-1].origin[1]:
                    tilted_chord_TE[i][1] = -(K[i][1] - K_TE[sect_index-1][1]) / (K[i][0] - K_TE[sect_index-1][0])
                    tilted_chord_TE[i][0] = K[i][0] * tilted_chord_TE[i][1] + K[i][1]
                  
            else:
                tilted_chord_TE[i][1] = -(K[i][1] - K_TE[sect_index][1]) / (K[i][0] - K_TE[sect_index][0])
                tilted_chord_TE[i][0] = K[i][0] * tilted_chord_TE[i][1] + K[i][1]

                if tilted_chord_TE[i][1] > wing.Segments[sect_index].origin[1]:
                    tilted_chord_TE[i][1] = -(K[i][1] - K_TE[sect_index+1][1]) / (K[i][0] - K_TE[sect_index+1][0])
                    tilted_chord_TE[i][0] = K[i][0] * tilted_chord_TE[i][1] + K[i][1]            
                
                tilted_chord_LE[i][1] = -(K[i][1] - K_LE[sect_index][1]) / (K[i][0] - K_LE[sect_index][0])
                tilted_chord_LE[i][0] = K[i][0] * tilted_chord_LE[i][1] + K[i][1]
                if tilted_chord_LE[i][1] < wing.Segments[sect_index-1].origin[1]:
                    tilted_chord_LE[i][1] = -(K[i][1] - K_LE[sect_index-1][1]) / (K[i][0] - K_LE[sect_index-1][0])
                    tilted_chord_LE[i][0] = K[i][0] * tilted_chord_LE[i][1] + K[i][1]

        tilted_chord[i] = ((tilted_chord_TE[i][0]-tilted_chord_LE[i][0])**2+(tilted_chord_TE[i][1]-tilted_chord_LE[i][1])**2)**0.5

    # Interpolate segments to find the strips airfoil
    num_airfoil_points = np.shape(wing.Segments[0].Airfoil.airfoil.points)
    Strip_points       = strip_airfoil_coordinates(wing,K,drag_sect,n_segments,num_airfoil_points[0],tilted_chord_TE,Y_strip)

    # Plot results
    plot_wing(num_airfoil_points,n_segments,wing,Strip_points,drag_sect,tilted_chord_LE,tilted_chord_TE,K_Trans,X_trans,Y_trans)

    return Strip_points, Y_strip, Trans_sweep, c, tilted_chord, drag_sect

def strip_airfoil_coordinates(wing,K,drag_sect,n_segments,num_airfoil_points,tilted_chord_TE,Y_strip):
    """Inerpolates airfoil strip coordinates along the wing 

    Assumptions:

    Source:

    Inputs:
    wing
    K
    drag_sect                                    [Unitless]
    n_segments                                   [Unitless]
    num_airfoil_points                           [Unitless]
    tilted_chord_TE                              [m]
    Y_strip                                      [m]
      
    Outputs:
    Strip_points                                 [m]


    Properties Used:
    N/A
    
    """
    
    # Find equations of the line for each airfoil point
    Ksect  = np.zeros((n_segments-1,num_airfoil_points,2))
    Strip_points = np.zeros((3,drag_sect,num_airfoil_points))
    Xstrip = np.zeros((drag_sect,num_airfoil_points))
    Ystrip = np.zeros((drag_sect,num_airfoil_points))
    Zstrip = np.zeros((drag_sect,num_airfoil_points))
    for i in range(n_segments-1):
        for j in range(num_airfoil_points):
            Ksect[i,j,0] = (wing.Segments[i+1].Airfoil.airfoil.points[j,0] - wing.Segments[i].Airfoil.airfoil.points[j,0])/  \
                           (wing.Segments[i+1].Airfoil.airfoil.points[j,1] - wing.Segments[i].Airfoil.airfoil.points[j,1])     
            Ksect[i,j,1] = wing.Segments[i+1].Airfoil.airfoil.points[j,0]  - Ksect[i,j,0]*wing.Segments[i+1].Airfoil.airfoil.points[j,1]
    
    # Find intersection points of strips and airfoil point lines
    for i in range(drag_sect):
        n = 1
        while tilted_chord_TE[i][1] > wing.Segments[n].Airfoil.airfoil.points[0,1]:
            n = n + 1
        for j in range(num_airfoil_points):
            # Check the special case - rectangular wing with zero sweep
            if K[i+1,0] == float('+inf') or K[i+1,0] == float('-inf'):
                Strip_points[1,i,j] = Y_strip[i+1]
                Strip_points[0,i,j] = Ksect[n-1,j,1]
                t = 0
                m = n
            else:   
                Strip_points[1,i,j] = (K[i+1,1] - Ksect[n-1,j,1]) / (Ksect[n-1,j,0]-K[i+1,0])
                Strip_points[0,i,j] = K[i+1,0] * Strip_points[1,i,j] + K[i+1,1]
                m = n
                if Strip_points[1,i,j] < wing.Segments[n-1].Airfoil.airfoil.points[0,1]:
                    Strip_points[1,i,j] = (K[i+1,1] - Ksect[n-2,j,1]) / (Ksect[n-2,j,0]-K[i+1,0])
                    Strip_points[0,i,j] = K[i+1,0] * Strip_points[1,i,j] + K[i+1,1]
                    m = n-1
                elif Strip_points[1,i,j] > wing.Segments[n].Airfoil.airfoil.points[0,1]:
                    Strip_points[1,i,j] = (K[i+1,1] - Ksect[n,j,1]) / (Ksect[n,j,0]-K[i+1,0])
                    Strip_points[0,i,j] = K[i+1,0] * Strip_points[1,i,j] + K[i+1,1]
                    m = n+1
                t = (Strip_points[0,i,j]-wing.Segments[m-1].Airfoil.airfoil.points[j,0])/  \
                (wing.Segments[m].Airfoil.airfoil.points[j,0]-wing.Segments[m-1].Airfoil.airfoil.points[j,0])

            Strip_points[2,i,j] = (wing.Segments[m].Airfoil.airfoil.points[j,2] - wing.Segments[m-1].Airfoil.airfoil.points[j,2])*t +  \
                           wing.Segments[m-1].Airfoil.airfoil.points[j,2]

    return Strip_points
        
def compute_strip_sweep(wing,Trans_sweep,Lift_distribution):  
    """Modifies initial conditions to the ones normal to the sweep line

    Assumptions:

    Source:

    Inputs:
    wing
    Trans_sweep                                  [Unitless]
    Lift_distribution                            Data structure with strip section properties
    
      
    Outputs:
    Lift_distribution                            Data structure with strip section properties
    
    Properties Used:
    N/A
    
    """

    n_cases = len(Lift_distribution.keys())
    for i in range(n_cases):
        Mach = Lift_distribution[i].Mach
        Lift_distribution[i].sweep = []
        for j in range(len(Lift_distribution[i].Ystrip)):
            sweep_index = 1
            while Lift_distribution[i].Ystrip[j] > 0.5 * wing.Segments[sweep_index].percent_span_location * wing.spans.projected:
                sweep_index = sweep_index + 1

            Lift_distribution[i].sweep.append(Trans_sweep[sweep_index])
            
    return Lift_distribution
   
def plot_wing(num_airfoil_points,n_segments,wing,Strip_points,drag_sect,tilted_chord_LE,tilted_chord_TE,K_Trans,X_trans,Y_trans):
    """Plots the wing planform with the airfoil strips required for the calculations

    Assumptions:

    Source:

    Inputs:
    num_airfoil_points
    n_segments
    Strip_points
    drag_sect
    tilted_chord_LE                              [m]
    tilted_chord_TE                              [m]
    K_Trans                                      [Unitless]
    X_trans                                      [m]
    Y_trans                                      [m]  
    wing

    Properties Used:
    N/A
    
    """    

    Nupper = 1
    Nlower = 0
    for i in range(1,num_airfoil_points[0]):
        if wing.Segments[0].Airfoil.airfoil.points[i-1,0] >= wing.Segments[0].Airfoil.airfoil.points[i,0]:
            Nupper = Nupper+1
        else:
            Nlower = Nlower+1
       
    X_sectu = np.zeros((Nupper,n_segments))
    Y_sectu = np.zeros((Nupper,n_segments))
    Z_sectu = np.zeros((Nupper,n_segments))
    X_sectl = np.zeros((Nlower,n_segments))
    Y_sectl = np.zeros((Nlower,n_segments))
    Z_sectl = np.zeros((Nlower,n_segments))
    for i in range(n_segments):
        n = 1
        m = 0
        X_sectu[0,i] = wing.Segments[i].Airfoil.airfoil.points[0,0]
        Y_sectu[0,i] = wing.Segments[i].Airfoil.airfoil.points[0,1]
        Z_sectu[0,i] = wing.Segments[i].Airfoil.airfoil.points[0,2]      
        for j in range(1,num_airfoil_points[0]):
            if wing.Segments[i].Airfoil.airfoil.points[j,0] <= wing.Segments[i].Airfoil.airfoil.points[j-1,0]:
                X_sectu[n,i] = wing.Segments[i].Airfoil.airfoil.points[j,0]
                Y_sectu[n,i] = wing.Segments[i].Airfoil.airfoil.points[j,1]
                Z_sectu[n,i] = wing.Segments[i].Airfoil.airfoil.points[j,2]
                n=n+1
            else:
                X_sectl[m,i] = wing.Segments[i].Airfoil.airfoil.points[j,0]
                Y_sectl[m,i] = wing.Segments[i].Airfoil.airfoil.points[j,1]
                Z_sectl[m,i] = wing.Segments[i].Airfoil.airfoil.points[j,2]
                m=m+1
                
   # plot the wing with chord sections
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(azim=0, elev=90)
    #fig, ax = plt.subplots()
    '''for i in range(1,n_segments):
        ax.plot([wing.Segments[i-1].origin[0],wing.Segments[i].origin[0]], \
                [wing.Segments[i-1].origin[1],wing.Segments[i].origin[1]], 'k')
        ax.plot([wing.Segments[i-1].Airfoil.airfoil.points[0][0],wing.Segments[i].Airfoil.airfoil.points[0][0]], \
                [wing.Segments[i-1].Airfoil.airfoil.points[0][1],wing.Segments[i].Airfoil.airfoil.points[0][1]], 'k')
        ax.plot([wing.Segments[i-1].origin[0],wing.Segments[i-1].Airfoil.airfoil.points[0][0]],   \
                 [wing.Segments[i-1].origin[1],wing.Segments[i-1].Airfoil.airfoil.points[0][1]], 'k')
        ax.plot([wing.Segments[i].origin[0],wing.Segments[i].Airfoil.airfoil.points[0][0]],       \
                 [wing.Segments[i].origin[1],wing.Segments[i].Airfoil.airfoil.points[0][1]], 'k')
    for i in range(1,drag_sect+1):
        ax.plot([tilted_chord_LE[i][0],tilted_chord_TE[i][0]],[tilted_chord_LE[i][1],tilted_chord_TE[i][1]], 'r')

    X_trans = np.zeros(n_segments)
    Y_trans = np.zeros(n_segments)
    Y_trans[0] = wing.Segments[0].origin[1]
    X_trans[0] = (Y_trans[0] - K_Trans[1][1]) / K_Trans[1][0]    
    for i in range(1,n_segments):
        Y_trans[i] = wing.Segments[i].origin[1]
        X_trans[i] = (Y_trans[i] - K_Trans[i][1]) / K_Trans[i][0]
    for i in range(drag_sect):
        ax.scatter(X_strip[i,:],Y_strip[i,:],color='g')
    ax.scatter(tilted_chord_LE[1:,0],tilted_chord_LE[1:,1])
    ax.set(xlabel='X, m', ylabel='Y, m', title='Wing Planform')
    ax.plot(X_trans,Y_trans)'''
    
    ax.plot_wireframe(X_sectu,Y_sectu,Z_sectu,color = 'k')
    ax.plot_wireframe(X_sectl,Y_sectl,Z_sectl,color = 'k')
    ax.scatter(Strip_points[0,:,:],Strip_points[1,:,:],Strip_points[2,:,:], color = 'g')
    #ax.set_xlim3d(-5, 40)
    #ax.set_ylim3d(-5, 40)
    #ax.set_zlim3d(-5, 40)
    plt.savefig('planform.png')

    return
