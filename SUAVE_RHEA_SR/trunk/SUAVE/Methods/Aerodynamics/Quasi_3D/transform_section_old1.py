## @ingroup Methods-Aerodynamics-Qasi_3D
#scale_airfoils.py
# 
# Created:  Aug 2019, S. Karpuk
# Modified: 
#

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

from SUAVE.Core import Data
import numpy as np
import math

## @ingroup Methods-Aerodynamics-Qasi_3D
def transform_section(geometry_file,wing):
    """ This functions scales, translates, and rotates given airfoil sections
        according to the .avl file

    Assumptions:
        None
        
    Source:
        Drela, M. and Youngren, H., AVL, http://web.mit.edu/drela/Public/web/avl

    Inputs:
        None

    Outputs:
        results     

    Properties Used:
        N/A
    """ 

    avl_section_parameters = []
    # read parameters from .avl file (X,Y,Z,chord,incidence)
    num_lines = len(open(geometry_file).readlines())
    with open(geometry_file, 'r') as f:
        for line in f:
            if 'main_wing' in line:
                num = 0
                while '#-----' not in line and num != num_lines:
                    line = f.readline()
                    if '#Xle' in line:
                        num_line = f.readline()
                        string = np.array(num_line.split())
                        avl_section_parameters.append(string.astype(np.float))
                    num = num + 1
    
    # transform airfoil coordinates to actual wing geometric coordinates
    if len(wing.Segments.keys())>0:

        n_segments  = len(wing.Segments.keys())
        print(n_segments)
        i_segs = 0
        for segment in wing.Segments.values():
            points = segment.Airfoil.airfoil.points
            transformed_points = np.zeros((np.shape(points)[0],3))

            # scale the section
            transformed_points[:,0] = points[:,0] * avl_section_parameters[i_segs][3]
            transformed_points[:,2] = points[:,1] * avl_section_parameters[i_segs][3]
            # rotate the section
            transformed_points[:,0] =  transformed_points[:,0] * math.cos(math.radians(avl_section_parameters[i_segs][4])) + \
                                       transformed_points[:,2] * math.sin(math.radians(avl_section_parameters[i_segs][4]))
            transformed_points[:,2] = -transformed_points[:,0] * math.sin(math.radians(avl_section_parameters[i_segs][4])) + \
                                       transformed_points[:,2] * math.cos(math.radians(avl_section_parameters[i_segs][4]))
            # translate the section
            transformed_points[:,0] = transformed_points[:,0] + avl_section_parameters[i_segs][0]
            transformed_points[:,1] = avl_section_parameters[i_segs][1]
            transformed_points[:,2] = transformed_points[:,2] + avl_section_parameters[i_segs][2]
            segment.Airfoil.airfoil.points = transformed_points
            segment.origin                 = transformed_points[np.argmin(transformed_points[:,0]),:]
            #print('here')
            print(i_segs)
            print(segment.tag)
            #print(wing.Segments[i_segs].Airfoil.airfoil.points)
            i_segs = i_segs+1
        #print(wing)
    
    return        


            
