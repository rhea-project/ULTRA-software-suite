#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
import numpy as np
import csv

import sys
#sys.path.append('C:/Uni/SE2A/system_methodology/03_tool_box/SUAVE/SUAVE-2.0.0/trunk')
#sys.path.append('D:/IfES/06_SE2A/system_methodology/03_tool_box/hensup')

# local application/library specific imports
#import SUAVE
#from SUAVE.Core.Data import Data
from SUAVE.OpenModelica.core.Data import Data

#-----------------------------------------------------------------------------#
# Superclass for atmospheric models
#-----------------------------------------------------------------------------#
class Profile(Data):
    def __init__(self):
        pass
        
    def getFileData(self, file_path=None):      
        if file_path is None:
            raise ValueError()
                    
        with open(file_path) as csvfile:
            reader = csv.DictReader(csvfile, delimiter=';')
            
            for field in reader.fieldnames:
                self[field] = []
            
            for row in reader:
                for field in reader.fieldnames:
                    self[field].append(float(row[field]))
            
        for field in reader.fieldnames:
            self[field] = np.asarray(self[field])
            
            