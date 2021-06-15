
from cases.hangar.richards_wing import Baseline
import numpy as np


def build_wing_data(sweep):

     u_inf = 28
     M = 4  # chordwise panels
     N = 11  # spanwise panels
     Msf = 5  # wake length in chord numbers
     rho_fact = 1.  # air density factor

     use_euler = True  # use euler angles or quaternions as orientation parametrisation
     if use_euler:
         orient = 'euler'
     else:
         orient = 'quat'
         
     case_notes = '_rho%s' % (str(rho_fact))
     case_rmks = 'M%gN%gMsf%g' % (M, N, Msf)

     # M4N11Msf5 - trim values
     alpha_deg = 4.5135
     cs_deflection = 0.1814
     thrust = 5.5129

     # M8N11Msf5 - trim values
     # alpha_deg = 4.5162
     # cs_deflection = 0.2373
     # thrust = 5.5129
     ws = Baseline(M=M,
              N=N,
              Mstarfactor=Msf,
              u_inf=u_inf,
              rho=1.02,
              alpha_deg=alpha_deg,
              roll_deg=0,
              cs_deflection_deg=cs_deflection,
              thrust=thrust,
              physical_time=20,
              case_name='horten_%s' % orient,
              case_name_format=4,
              case_remarks=case_rmks + case_notes)

     ws.set_properties()
     ws.sweep_LE = sweep
     ws.initialise()
     ws.clean_test_files()

     # ws.update_mass_stiffness(sigma=1., sigma_mass=1.5)
     ws.update_mass_stiffness(sigma=1., sigma_mass=2.5)
     ws.update_fem_prop()
     ws.generate_fem_file()
     ws.update_aero_properties()
     ws.generate_aero_file()

     return ws
