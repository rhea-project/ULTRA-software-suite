## @ingroup Methods-Aerodynamics-Xfoil
#create_avl_datastructure.py
# 
# Created:  Jul 2019, S. Karpuk
# Modified: 
#    


# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------
import scipy
import numpy as np
import math
import matplotlib.pyplot as plt

from copy import deepcopy

# SUAVE Imports
import SUAVE
from .Data.Inputs   import Inputs

# SUAVE-Xfoil Imports
from SUAVE.Core import Data, Units

## @ingroup Methods-Aerodynamics-Xfoil
def create_xfoil_datastructure(geometry,ref_line,wing):
        """ This translates the aircraft geometry into the format used in the AVL run file

        Assumptions:
            None
    
        Source:
            Drela, M. and Youngren, H., AVL, http://web.mit.edu/drela/Public/web/avle
    
        Inputs:
            geometry    
    
        Outputs:
            avl_inputs
    
        Properties Used:
            N/A
        """    
        xfoil_airfoils             = translate_xfoil_geometry(geometry,ref_line,wing)

        # pack results in a new Xfoil inputs structure
        xfoil_inputs               = Inputs()
        xfoil_inputs.airfoils      = xfoil_airfoils

        return xfoil_inputs


def translate_xfoil_geometry(geometry,ref_line,wing):
        """ Translates geometry from the vehicle setup to Xfoil format

        Assumptions:
            None

        Source:
            None

        Inputs:
            geometry
                geometry.airfoils                                                [data stucture] 


        Outputs:
            airfoils - aircraft geometry in AVL format                           [data stucture] 

        Properties Used:
            N/A
        """ 

        Xfoil_airfoils = SUAVE.Components.Wings.Main_Wing()
        airfoil_data_size = np.shape(geometry)

        # Create Data structure of each airfoil
        incidence = np.zeros(airfoil_data_size[1])
        for i in range(airfoil_data_size[1]):
                segment = SUAVE.Components.Wings.Segment()
                strip   = SUAVE.Components.Wings.Airfoils.Airfoil()
                # Find the airfoil incidence and chord length
                Xmax           = max(geometry[0,i,:])
                Xmax_ind       = np.argmax(geometry[0,i,:])
                Xmin           = min(geometry[0,i,:])
                Xmin_ind       = np.argmin(geometry[0,i,:])
                LE             = [Xmin,geometry[1,i,Xmin_ind],geometry[2,i,Xmin_ind]]
                TE             = [Xmax,geometry[1,i,Xmax_ind],geometry[2,i,Xmax_ind]]
                chord          = ((TE[0]-LE[0])**2+(TE[1]-LE[1])**2+(TE[2]-LE[2])**2)**0.5
                incidence[i]   = -math.atan((TE[2]-LE[2])/((TE[0]-LE[0])**2+(TE[1]-LE[1])**2)**0.5)
                #print(LE,TE,((TE[0]-LE[0])**2+(TE[1]-LE[1])**2)**0.5,TE[2]-LE[2],incidence[i])
                ref_line_point = [(TE[0]-LE[0])*ref_line+LE[0], (TE[1]-LE[1])*ref_line+LE[1], (TE[2]-LE[2])*ref_line+LE[2]]

                strip_points = np.zeros((airfoil_data_size[2],3))
                points       = np.zeros((airfoil_data_size[2],2))

                #Transform the strip airfoil into a 2D normilized format
                sect = 0
                while ref_line_point[1] > wing.Segments[sect].origin[1] :
                        sect = sect + 1
                        
                sweep = wing.Segments[sect-1].sweeps.quarter_chord
                for j in range(airfoil_data_size[2]):                    
                        # Translate the coordinate axis to the local strip LE
                        strip_points[j,:] = [geometry[0,i,j]-LE[0], geometry[1,i,j]-LE[1], geometry[2,i,j]-LE[2]]
                       
                        points[j,0] = math.cos(sweep)*math.cos(incidence[i])*strip_points[j,0] - \
                                      math.sin(sweep)*math.cos(incidence[i])*strip_points[j,1] - \
                                      math.sin(incidence[i])*strip_points[j,2]
                        test_point  = math.sin(sweep)*strip_points[j,0] + math.cos(sweep)*strip_points[j,1]
                        points[j,1] = math.cos(sweep)*math.sin(incidence[i])*strip_points[j,0] - \
                                      math.sin(sweep)*math.sin(incidence[i])*strip_points[j,1] + \
                                      math.cos(incidence[i])*strip_points[j,2]

                # Normalize airfoils
                points = np.divide(points, chord)
                
                # plot airfoils
                '''fig, ax = plt.subplots()
                ax.scatter(points[:,0],points[:,1])
                ax.set_xlim(0,1)
                ax.set_ylim(-0.5,0.5)
                plt.show()'''

                # Fill the Data structure
                segment.tag           = 'Airfoil_' + str(i)
                segment.chord         = chord
                segment.twist         = incidence[i]
                strip.coordiante_file = segment.tag + '.dat'
                strip.points          = points
                segment.append_airfoil(strip)
                Xfoil_airfoils.Segments.append(segment)

        return Xfoil_airfoils



