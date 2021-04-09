## @ingroup Methods-Weights-EMWET
# main_wing.py
#
# Created:  Jun 2020, S. Karpuk
# Modified:   

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import sys
import time
import subprocess
import SUAVE
from SUAVE.Core import Units, Data, redirect
from SUAVE.Analyses.Atmospheric                                       import US_Standard_1976 as atmosphere
from SUAVE.Methods.Aerodynamics                                       import AVL_Loads
from SUAVE.Methods.Aerodynamics.AVL.purge_files                       import purge_files

from SUAVE.Methods.Geometry.Three_Dimensional import Compute_sweep
import numpy as np
import scipy.interpolate as interpol
import os
import re

# ----------------------------------------------------------------------
#   Wing weight using EMWET
# ----------------------------------------------------------------------

## @ingroup Methods-Weights-EMWET-Tube_Wing_EMWET
def wing_main(vehicle,eng_origin,num_eng_w,wt_engine):      
    """ Calculate the weight of the wing using the EMWET method
    
    Assumptions:
        Now wing twist is considered for the wingbox simulation
    
    Source: 
        A. Elham, G. La Rocca, M.J.L van Tooren, Development and implementation of an advanced,
        design-sensitive method for wing weight estimation, Aerospace Science and Technology, 29 (2013), 100-113
        
    Inputs:
        vehicle.num_eng_w                      [Unitless]
               .eng_origin                     [Unitless]

    Outputs:
        weight - weight of the horizontal tail                                                 [kilograms]
       
    Properties Used:
        N/A
    """

    # Set constants up
    call_EMWET = r'C:\Users\May\Desktop\SUAVE_TUBS_se2a_V1\trunk\SUAVE\Methods\Weights\EMWET\EMWET.exe'
    log_file   = r'EMWET.log'
    err_file   = r'EMWET.err'

    # Create the initiation file
    write_init_file(vehicle,num_eng_w,eng_origin,wt_engine)

    # Create a loads file
    limit_load_distribution(vehicle)

    # Run EMWET
    Run_EMWET(vehicle._base.tag,call_EMWET,log_file,err_file)

    # Read the EMWET output weight
    weight = read_weight(vehicle)


    return weight


def write_init_file(vehicle,num_eng_w,eng_origin,wt_engine):
    """ Creates an initiation fole for EMWET
    
    Assumptions:
   
    Source: 
    
    Inputs:
        wing_Segments
            percent_span_location               [dimensionless]
            root_chord_percent                  [dimensionless]
            dihedral_outboard                   [dimensionless]
            sweeps
                quarter_chord                   [dimensionless]
            front_spar                          [dimensionless]
            rear_spar                           [dimensionless]
            airfoil        
            
        E          - an array of Young's muduli                                                [N/meters**2]
        sigma_maxt - ultimate tensile strength                                                 [N/meters**2]
        sigma_maxc - ultimate compressive strength                                             [N/meters**2]
        rho        - material density                                                          [kg/meters**3]
        rib_pitch  - rib pitch                                                                 [meters]
        eta_pan    - panel efficiency                                                          [dimensionless]
        Nult       - ultimate design load of the aircraft                                      [dimensionless]
        TOW        - maximum takeoff weight of the aircraft                                    [kilograms]
        wt_zf      - maximum zero-fuel weight                                                  [kilograms]
        wt_prop    - propulsion system weight                                                  [kilograms]
        eng_origin - an array of engine origins                                                [meters]
        S_ref      - wing reference area                                                       [meters**2]
        b          - wing span                                                                 [meters]
        fuel tank  - array of fuel tank reference locations                                    [dimensionless]

    Outputs:
        weight - weight of the horizontal tail                                                 [kilograms]
       
    Properties Used:
        N/A
    """

    # Unpack inputs
    TOW  = vehicle.mass_properties.max_takeoff
    ZFW  = vehicle.mass_properties.max_zero_fuel  
    Nlim = vehicle.envelope.limit_load 
    Sref = vehicle.wings['main_wing'].areas.reference
    b    = 0.5 * vehicle.wings['main_wing'].spans.projected
    root = vehicle.wings['main_wing'].chords.root
    wing_Segments = vehicle.wings['main_wing'].Segments

    n_segments = len(wing_Segments.keys())
    geom_perc_span          = np.zeros(n_segments)
    geom_root_chord_percent = np.zeros(n_segments)
    geom_dihedral_outboard  = np.zeros(n_segments)
    geom_sweep_quart        = np.zeros(n_segments)
    geom_front_spar         = np.zeros(n_segments)
    geom_rear_spar          = np.zeros(n_segments)
    airfoil                 = []
    
    for i_segm in range(n_segments):
        geom_perc_span[i_segm]          = wing_Segments[i_segm].percent_span_location 
        geom_root_chord_percent[i_segm] = wing_Segments[i_segm].root_chord_percent 
        geom_dihedral_outboard[i_segm]  = wing_Segments[i_segm].dihedral_outboard
        geom_sweep_quart[i_segm]        = wing_Segments[i_segm].sweeps.quarter_chord
        geom_front_spar[i_segm]         = wing_Segments[i_segm].front_spar  
        geom_rear_spar[i_segm]          = wing_Segments[i_segm].rear_spar
        airfoil.append(wing_Segments[i_segm].Airfoil.airfoil.coordinate_file)

    # Unpack material properties
    E   = np.zeros(4)
    rho = np.zeros(4)
    material_flag = vehicle.wings['main_wing'].materials.type
    if material_flag == 'composite':
        m = np.zeros(4)
        n = np.zeros(4)

        E[0] = vehicle.wings['main_wing'].materials.upper_panel.E
        E[1] = vehicle.wings['main_wing'].materials.lower_panel.E
        E[2] = vehicle.wings['main_wing'].materials.front_spar.E
        E[3] = vehicle.wings['main_wing'].materials.rear_spar.E

        rho[0] = vehicle.wings['main_wing'].materials.upper_panel.density
        rho[1] = vehicle.wings['main_wing'].materials.lower_panel.density
        rho[2] = vehicle.wings['main_wing'].materials.front_spar.density
        rho[3] = vehicle.wings['main_wing'].materials.rear_spar.density   

        m[0] = vehicle.wings['main_wing'].materials.upper_panel.m
        m[1] = vehicle.wings['main_wing'].materials.lower_panel.m
        m[2] = vehicle.wings['main_wing'].materials.front_spar.m
        m[3] = vehicle.wings['main_wing'].materials.rear_spar.m

        n[0] = vehicle.wings['main_wing'].materials.upper_panel.n
        n[1] = vehicle.wings['main_wing'].materials.lower_panel.n
        n[2] = vehicle.wings['main_wing'].materials.front_spar.n
        n[3] = vehicle.wings['main_wing'].materials.rear_spar.n

        
    else:   
        sigma_skin = np.zeros((2,2))
        sigma_spar = np.zeros([7,2])
                
        E[0] = vehicle.wings['main_wing'].materials.upper_panel.E
        E[1] = vehicle.wings['main_wing'].materials.lower_panel.E
        E[2] = vehicle.wings['main_wing'].materials.front_spar.E
        E[3] = vehicle.wings['main_wing'].materials.rear_spar.E      
        
        sigma_skin[0,0] = vehicle.wings['main_wing'].materials.upper_panel.sigma_maxt
        sigma_skin[0,1] = vehicle.wings['main_wing'].materials.lower_panel.sigma_maxt
        sigma_skin[1,0] = vehicle.wings['main_wing'].materials.upper_panel.sigma_maxc
        sigma_skin[1,1] = vehicle.wings['main_wing'].materials.lower_panel.sigma_maxc

        sigma_spar[0,0] = vehicle.wings['main_wing'].materials.front_spar.longitudinal.sigma_maxt
        sigma_spar[1,0] = vehicle.wings['main_wing'].materials.front_spar.longitudinal.sigma_maxc
        sigma_spar[2,0] = vehicle.wings['main_wing'].materials.front_spar.longitudinal.sigma_ultt
        sigma_spar[3,0] = vehicle.wings['main_wing'].materials.front_spar.lateral.sigma_maxt
        sigma_spar[4,0] = vehicle.wings['main_wing'].materials.front_spar.lateral.sigma_maxc
        sigma_spar[5,0] = vehicle.wings['main_wing'].materials.front_spar.lateral.sigma_ultt
        sigma_spar[6,0] = vehicle.wings['main_wing'].materials.front_spar.sigma_maxs

        sigma_spar[0,1] = vehicle.wings['main_wing'].materials.rear_spar.longitudinal.sigma_maxt
        sigma_spar[1,1] = vehicle.wings['main_wing'].materials.rear_spar.longitudinal.sigma_maxc
        sigma_spar[2,1] = vehicle.wings['main_wing'].materials.rear_spar.longitudinal.sigma_ultt
        sigma_spar[3,1] = vehicle.wings['main_wing'].materials.rear_spar.lateral.sigma_maxt
        sigma_spar[4,1] = vehicle.wings['main_wing'].materials.rear_spar.lateral.sigma_maxc
        sigma_spar[5,1] = vehicle.wings['main_wing'].materials.rear_spar.lateral.sigma_ultt
        sigma_spar[6,1] = vehicle.wings['main_wing'].materials.rear_spar.sigma_maxs
                
        rho[0] = vehicle.wings['main_wing'].materials.upper_panel.density
        rho[1] = vehicle.wings['main_wing'].materials.lower_panel.density
        rho[2] = vehicle.wings['main_wing'].materials.front_spar.density
        rho[3] = vehicle.wings['main_wing'].materials.rear_spar.density
           
    rib_pitch  = vehicle.wings['main_wing'].structures.rib_pitch
    eta_pan    = vehicle.wings['main_wing'].structures.panel_efficiency 

    # Unpack fuel tank properties
    fuel_tank    = [vehicle.wings['main_wing'].fuel_tank.start, vehicle.wings['main_wing'].fuel_tank.end]  

    # Unpack engines
    wing_eng = []
    wing_eng_num = 0
    eng_size = np.shape(eng_origin)[0]
    for i in range(eng_size):
        if eng_origin[i][1]>0:
            wing_eng_num += 1
            wing_eng.append(eng_origin[i][1])

    # Calculate wing characteristics
    XLE    = np.zeros(n_segments)
    YLE    = np.zeros(n_segments) 
    ZLE    = np.zeros(n_segments)  
    chords = np.zeros(n_segments)  

    # Shift the reference to the wing root LE
    XLE[0] = 0;     YLE[0] = 0;     ZLE[0] = 0;   chords[0] = root  
    for i in range(1,n_segments):
        YLE[i]    = geom_perc_span[i] * b
        chords[i] = geom_root_chord_percent[i] * root

        dY        = (geom_perc_span[i]-geom_perc_span[i-1]) * b
        LE_sweep  = Compute_sweep.Compute_sweep(dY,chords[i-1],chords[i],0.25,geom_sweep_quart[i-1],'to leading')

        XLE[i] = XLE[i-1] + dY*np.tan(LE_sweep)
        ZLE[i] = ZLE[i-1] + dY*np.tan(geom_dihedral_outboard[i-1])

    #write input into a file
    f = open(vehicle._base.tag + ".init", "w")
    f.write(str(TOW) + ' ' + str(ZFW) + '\n')
    f.write(str(Nlim) + '\n')
    f.write(str(Sref) + ' ' + str(2*b) + ' ' + str(n_segments) + ' ' + str(n_segments) + '\n')
    for i in range(n_segments):
        f.write(str(geom_perc_span[i]) + ' ' + str(airfoil[i]) + '\n')
    for i in range(n_segments):
        f.write(str(chords[i]) + ' ' + str(XLE[i]) + ' ' + str(YLE[i]) + ' ' + str(ZLE[i]) + ' ' +   \
            str(geom_front_spar[i]) + ' ' + str(geom_rear_spar[i]) + '\n')
    f.write(str(fuel_tank[0])  + ' ' + str(fuel_tank[1])+ '\n')
    f.write(str(0.5*wing_eng_num) + '\n')         # Assume each wing side has an equal amount of engines
    for i in range(int(0.5*wing_eng_num)):
        f.write(str(wing_eng[2*i]/b) + ' ' + str(wt_engine[2*i]) + '\n')
    
    f.write(material_flag + '\n')
    if material_flag == 'composite':
        f.write(str(E[0]) + ' ' +str(rho[0]) + ' ' +str(m[0]) + ' ' +str(n[0]) + '\n')
        f.write(str(E[1]) + ' ' +str(rho[1]) + ' ' +str(m[1]) + ' ' +str(n[1]) + '\n')
        f.write(str(E[2]) + ' ' +str(rho[2]) + ' ' +str(m[2]) + ' ' +str(n[2]) + '\n')
        f.write(str(E[3]) + ' ' +str(rho[3]) + ' ' +str(m[3]) + ' ' +str(n[3]) + '\n')
    else:
        f.write(str(E[0]) + ' ' +str(rho[0]) + ' ' +str(sigma_skin[0,0]) + ' ' +str(sigma_skin[1,0]) + '\n')
        f.write(str(E[1]) + ' ' +str(rho[1]) + ' ' +str(sigma_skin[0,1]) + ' ' +str(sigma_skin[1,1]) + '\n')
        f.write(str(E[2]) + ' ' +str(rho[2]) + ' ' +str(sigma_spar[0,0]) + ' ' +str(sigma_spar[1,0]) + ' ' +str(sigma_spar[3,0]) + ' ' + \
                str(sigma_spar[4,0]) + ' ' +str(sigma_spar[2,0]) + ' ' +str(sigma_spar[5,0]) + ' ' +str(sigma_spar[6,0]) + '\n')
        f.write(str(E[3]) + ' ' +str(rho[3]) + ' ' +str(sigma_spar[0,1]) + ' ' +str(sigma_spar[1,1]) + ' ' +str(sigma_spar[3,1]) + ' ' + \
                str(sigma_spar[4,1]) + ' ' +str(sigma_spar[2,1]) + ' ' +str(sigma_spar[5,1]) + ' ' +str(sigma_spar[6,1]) + '\n')
    f.write(str(eta_pan) + ' ' + str(rib_pitch)  + '\n')
    f.write(str(0))

    f.close()

    return

def limit_load_distribution(vehicle):
    """Runs AVL and obtains lift distribution under the limit load conditions

    Assumptions:

    Source:

    Inputs:
    vehicle
      
    Outputs:

    Properties Used:
    N/A
    
    """

    #   Initialize Cruise Configurations for the limit load case
    # ------------------------------------------------------------------
    TOW  = vehicle.mass_properties.max_takeoff
    Nlim = vehicle.envelope.limit_load 
    Sref = vehicle.wings['main_wing'].areas.reference
    b    = vehicle.wings['main_wing'].spans.projected
    wing_Segments = vehicle.wings['main_wing'].Segments
    n_segments    = len(wing_Segments.keys())
    sweep         = vehicle.wings['main_wing'].sweeps.quarter_chord
  
    geom_span = np.zeros(n_segments)
    for i_segm in range(n_segments):     
        geom_span[i_segm] = 0.5 * wing_Segments[i_segm].percent_span_location * b

    # Set the max load case flight conditions up 
    #---------------------------------------------------------------------------------------------
    atmo      = atmosphere()
    atmo_data = atmo.compute_values(0)

    avl                = Data()
    avl.current_status = Data()
    avl.current_status.batch_index = 0
    
    conditions = Data()
    conditions.weights      = Data()
    conditions.freestream   = Data()
    conditions.aerodynamics = Data()
    conditions.weights.total_mass     = TOW
    conditions.freestream.density     = atmo_data.density[0][0]
    conditions.freestream.gravity     = 9.81


    # Estimate CLmax for a cruise configuration (a simplified method without viscous corrections)
    #--cl max based on airfoil t_c
    Cl_max_ref = 0.9*1.5*np.cos(sweep)# -0.0009*tc**3 + 0.0217*tc**2 - 0.0442*tc + 0.7005

    w_Clmax = Cl_max_ref 

    V = (2*Nlim*TOW*9.81/(atmo_data.density[0][0]*Sref*w_Clmax))**0.5
    q = 0.5*atmo_data.density[0][0]*(V)**2 
    

    conditions.aerodynamics.lift_coefficient    = np.zeros(1)
    conditions.aerodynamics.lift_coefficient[0] = w_Clmax
    conditions.freestream.mach_number = V/atmo_data.speed_of_sound[0][0]

    # Import AVL settings
    avl.settings = AVL_Loads.Data.Settings()

    # Set spanwise and chordwise panel distribution
    avl.settings.discretization.defaults.wing.spanwise_vortices  = 25
    avl.settings.discretization.defaults.wing.chordwise_vortices = 10  
    
    # Add the geometry to AVL
    avl.geometry = vehicle
    
    # Assign settings
    run_folder                       = os.path.abspath(avl.settings.filenames.run_folder)
    output_template                  = avl.settings.filenames.output_template
    Trefftz_template                 = avl.settings.filenames.Trefftz_output_template
    batch_template                   = avl.settings.filenames.batch_template
    deck_template                    = avl.settings.filenames.deck_template

    # rename default avl aircraft tag
    avl.settings.filenames.features = vehicle.tag + '.avl'

    batch_index                     = avl.current_status.batch_index
    avl.current_status.batch_file   = batch_template.format(batch_index)
    avl.current_status.deck_file    = deck_template.format(batch_index)

    # translate conditions
    cases                           = AVL_Loads.translate_data.translate_conditions_to_cases(avl,conditions)
    avl.current_status.cases        = cases

    # Set the regression flag to False
    avl.regression_flag                 = False
    
    # case filenames
    for case in cases:
        cases[case].result_filename       = output_template.format(case)
        cases[case].loads_result_filename = Trefftz_template.format(case)

    # write the input files
    with redirect.folder(run_folder,force=False):
        AVL_Loads.write_geometry(avl)
        AVL_Loads.write_run_cases(avl)
        AVL_Loads.write_input_deck(avl)
    
        # RUN AVL!
        results_avl = AVL_Loads.run_analysis(avl)

    # open the loads file and read shear and moment distributions
    Y   = np.zeros(avl.settings.discretization.defaults.wing.spanwise_vortices)
    c   = np.zeros(avl.settings.discretization.defaults.wing.spanwise_vortices)
    Vz  = np.zeros(avl.settings.discretization.defaults.wing.spanwise_vortices)
    Mx  = np.zeros(avl.settings.discretization.defaults.wing.spanwise_vortices)

    lines_to_read = np.arange(start=20, stop=20+avl.settings.discretization.defaults.wing.spanwise_vortices, step=1)

    with redirect.folder(run_folder,force=False):
        for case in cases:
            f = open(cases[case].loads_result_filename)
            all_lines = f.readlines()
            for i in range(len(lines_to_read)):
                Y[i]   = 2*float(all_lines[lines_to_read[i]][8:16].strip()) / b
                c[i]   = float(all_lines[lines_to_read[i]][17:25].strip())
                Vz[i]  = float(all_lines[lines_to_read[i]][34:42].strip()) * q
                Mx[i]  = float(all_lines[lines_to_read[i]][89:97].strip()) * c[i]**2 * q
    # Correct for the root and tip of the wing
    Y[0] = 0;         Y[len(Y)-1] = 1

    # Write the loads file
    f = open(vehicle._base.tag + ".load", "w")
    for i in range(len(Y)):
        f.write(str(Y[i]) + ' ' + str(Vz[i]) + ' ' + str(Mx[i]) + '\n')
    f.close()

    return 

def Run_EMWET(case_tag,call_EMWET,log_file,err_file):
    """run EMWET 

    Assumptions:

    Source:

    Inputs:
        call_EMWET                 

    Outputs:

    Properties Used:
    N/A
    
    """

    if isinstance(log_file,str):
        purge_files(log_file)
    if isinstance(err_file,str):
        purge_files(err_file)
    
    with redirect.output(log_file,err_file):
    
        ctime = time.ctime() # Current date and time stamp

        print_output = False
                
        # Initialize suppression of console window output
        if print_output == False:
            devnull = open(os.devnull,'w')
            sys.stdout = devnull       
                    
        # Run AVL
        avl_run = subprocess.Popen([call_EMWET,case_tag],stdout=sys.stdout,stderr=sys.stderr,stdin=subprocess.PIPE)
        #for line in commands:
        #    avl_run.stdin.write(line.encode('utf-8'))
        #    avl_run.stdin.flush()
                  
        # Terminate suppression of console window output  
        if print_output == False:
            sys.stdout = sys.__stdout__                    
                    
        avl_run.wait()
    
        exit_status = avl_run.returncode
        ctime = time.ctime()

    return exit_status

def read_weight(vehicle):
    """read wing weight and wing section properties data

    Assumptions:
        Shear modulus is taken from the upper skin 
    Source:

    Inputs:
    vehicle
      
    Outputs:
    weight                                 [m]


    Properties Used:
    N/A
    
    """
    # Unpack inputs
    G   = vehicle.wings['main_wing'].materials.upper_panel.G
    TOW = vehicle.mass_properties.max_takeoff * 2.20462          # in lb

    count = len(open(vehicle._base.tag + '.weight').readlines(  ))
    Y     = np.zeros(count)
    J     = np.zeros(count)
    fspar = np.zeros(count)
    rspar = np.zeros(count)
    x0    = np.zeros(count)

    f   = open(vehicle._base.tag + '.weight', "r")
    all_lines = f.readlines()
    weight    = [float(s) for s in re.findall(r'-?\d+\.?\d*', all_lines[0])]
    CG        = [float(s) for s in re.findall(r'-?\d+\.?\d*', all_lines[2])]

    lines_to_read = np.arange(start=7, stop=count, step=1)
    for i in range(len(lines_to_read)):
        Y[i]     = float(all_lines[lines_to_read[i]][0:11].strip())  
        J[i]     = float(all_lines[lines_to_read[i]][95:108].strip())
        fspar[i] = float(all_lines[lines_to_read[i]][111:124].strip())
        rspar[i] = float(all_lines[lines_to_read[i]][127:140].strip())
        x0[i]    = float(all_lines[lines_to_read[i]][176:189].strip())

    for i in range(len(Y)-1):
        if Y[i] == Y[i-1]:
            Y[i]     = 0.5*(Y[i-1] + Y[i+1])
            J[i]     = 0.5*(J[i-1] + J[i+1])
            fspar[i] = 0.5*(fspar[i-1] + fspar[i+1])
            rspar[i] = 0.5*(fspar[i-1] + fspar[i+1])
            x0[i]    = 0.5*(x0[i-1] + x0[i+1])
    Jinterp     = interpol.interp1d(Y,J)
    x0interp    = interpol.interp1d(Y,x0)
    fsparinterp = interpol.interp1d(Y,fspar)
    rsparinterp = interpol.interp1d(Y,rspar)

    # Calculate section properties required for flutter analysis
    vehicle.wings['main_wing'].aeroelasticity.CG        = CG
    vehicle.wings['main_wing'].aeroelasticity.GJ        = G * J[0]
    vehicle.wings['main_wing'].aeroelasticity.GJ_mid    = G * Jinterp(0.5)
    vehicle.wings['main_wing'].aeroelasticity.x0_60     = x0interp(0.6)
    vehicle.wings['main_wing'].aeroelasticity.fspar_60  = fsparinterp(0.6)
    vehicle.wings['main_wing'].aeroelasticity.rspar_60  = rsparinterp(0.6)

    # Pack output
    vehicle.wings['main_wing'].MOI = Data()
    vehicle.wings['main_wing'].MOI.Polar_MOI = J

    return weight           # reduces composite wing