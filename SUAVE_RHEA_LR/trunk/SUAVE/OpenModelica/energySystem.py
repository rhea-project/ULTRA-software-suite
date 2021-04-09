#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports
import os

# related third party imports
import numpy as np
import scipy.io
from scipy.optimize import differential_evolution
from OMPython import ModelicaSystem

# local application/library specific imports
import SUAVE
from SUAVE.OpenModelica.de.dea import DifferentialEvolution
from SUAVE.OpenModelica.profile.flight_profile import Profile
from SUAVE.OpenModelica.core.Data import Data


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class EnergyModel(object):
    def __init__(self):
        self.DEF_PATH = hensup.__file__.replace('__init__.py','').replace("\\", "/")
        self.SRC_PATH = self.DEF_PATH + '/src/hensup'
        self.LIB_PATH = self.DEF_PATH + '/lib'
        self.DAT_PATH = self.DEF_PATH + '/data'
        self.MOD_PATH = self.DEF_PATH + '/input/models/aircraft_simulation'
        self.MOD_FILE = '/Aircraft_Turbofan.mo'
        self.MOD_NAME = 'Aircraft_Turbofan'
        self.RES_NAME = 'Aircraft_Turbofan_res.mo'
        
        self.MOD_FILE_PATH = self.MOD_PATH + self.MOD_NAME
        
        self.input_sim_vars = ['altitude', 'MachNumber', 'thrust']
        self.input_opt_vars = ['altitude', 'MachNumber', 'thrust']
        
        self.state = 0
        self.name = []
        
        self.delta = 100
        self.penalty = 1e+20
        
        self.inputs = Data()
        self.outputs = Data()
        
        self.profile = Profile()
        self.start_time = -1
        self.final_time = -1
        
        self.maximum_thrust = 0
                
        # inti optimisation attributes and objects
        self.objective = ['Aircraft_Model.massDesign',
                          'Aircraft_Model.deficit']

        self.designParameter = []
        
        self.designParameter_nb = len(self.designParameter)
        
        self.designParameter_bounds = []
        
        self.de = DifferentialEvolution(self.fobj, 
                                        self.designParameter_bounds)
                    
    def setTestExample(self):
        self.profile.getFileData(self.DAT_PATH + '/mission.csv')
        
        self.input_object = self.convertInput()
        
        self.start_time = self.input_object[self.input_sim_vars[0]][0][0]
        self.final_time = self.input_object[self.input_sim_vars[0]][-1][0]
        
        self.setInput()
        
        step_size = (self.final_time - self.start_time) / 100
            
        self.mod.setSimulationOptions(startTime=self.start_time, 
                                        stopTime=self.final_time, 
                                        stepSize=step_size, 
                                        tolerance=1e-06)
        
    def loadModel(self):
        self.mod = ModelicaSystem(self.MOD_FILE_PATH, self.MOD_NAME, ["lib"])
        
    def setDesignParameters(self, x):
        d = {}        
        for i in range(0,len(x)):
            d[self.designParameter[i]] = x[i]
        self.mod.setParameters(**d)
                
    def optimise(self, max_iter=10):
        self.de.its = max_iter
        res = list(self.de.de())
        
        return res        
    
    def convertInput(self):      
        # unpack
        inputs = self.inputs
        
        input_object = {}
        input_traj = {}
        
        inputs.time = self.profile.time
        for i in self.input_sim_vars:
            inputs[i] = self.profile[i]
            input_object[i] = []
            input_traj[i] = []
        
        for key in input_object:
            if inputs.time.shape[0] == 1:
                input_traj[key] = np.transpose(np.vstack((inputs.time, 
                          inputs[key])))
            else:
                input_traj[key] = np.vstack((np.transpose(inputs.time), 
                                             np.transpose(inputs[key])))
                
            for i in range(len(inputs.time)):
                input_object[key].append((input_traj[key][0][i], input_traj[key][1][i]))
        
        return input_object
    
    def setInput(self):
        self.mod.setInputs(**self.input_object)
        
    def simulate(self):
        return self.mod.simulate()
    
    def get_res(self, mat, pre, dataSet, pointer):
        return pre*mat[dataSet][pointer-1][:]
    
    def read_mat(self, filename, filepath=''):
        ''' 
        Auswertungsfunktion für *.mat-Files, die mit OMPython generiert wurden
        
        Eingangsvariablen:
            filename := Name des *.mat-File (String)
            filepath := Pfad zum *.mat-File (String) 
        Ausgangsvariablen:
            res := Dictonary, key := Variablenname des Modells (Parameter bekommen
                                                                nur zwei Werte 
                                                                übergeben)
        
        '''
        ###########################################################################
        # Einlesen der *.mat-Datei
        if len(filepath)>0 and filepath[-1] != '/':
            filename='/'+filename
        try:
            mat = scipy.io.loadmat(filepath+filename)
        except:
            print('This file does not exsist: ', filepath+filename)
            
        ###########################################################################
        # Auswerten der *.mat-Datei
        num_var = len(mat['name'][0])
        max_len = len(mat['name'])
        res = {}
        
        if self.state == 0:
            for k in range(num_var):
                i = 0
                # Einlesen der Variablennamen
                name = ''
                while i < max_len-1:
                    try:
                        if mat['name'][i][k] != '\x00':
                            name +=  mat['name'][i][k]
                            i += 1
                        else:
                            break
                    except:
                        break
                else:
                    pass
                self.name.append(name)
            self.state = 1
        else:
            pass
            
        for n in range(num_var):
            # Parameter:= data_1 / Variable := data_2
            if mat['dataInfo'][0][n] == 0:
                dataSet = 'data_2'
            else:
                dataSet = 'data_' + str(mat['dataInfo'][0][n])
                
            # Identifizieren von Alias
            pointer = mat['dataInfo'][1][n] 
            if pointer < 0:
                pointer *= -1
                pre = -1
            else:
                pre = 1
            res.update({self.name[n]: self.get_res(mat, pre, dataSet, pointer)})
        
        return res

    def fobj(self, x):
        
        self.setDesignParameters(x)
        
        try:
            self.simulate()
            
            res = self.read_mat(self.RES_NAME)
            
            mass_system = res[self.objective[0]][-1]
            thrustDeficit = res[self.objective[1]][-1]
            
            time = res['time']
            time_steps = np.diff(time)
            massFlowKerosene = res['Aircraft_Model.EMS_Kerosene.massFlowFuelSet'][-1]
            kerosene_mass = np.dot(time_steps, massFlowKerosene)
            
            mass = kerosene_mass + mass_system
            
            if thrustDeficit > self.delta:
                objective_value = self.penalty + mass + self.delta
            else:
                objective_value = thrustDeficit + mass
        except:
            objective_value = self.penalty
        
        return objective_value
    