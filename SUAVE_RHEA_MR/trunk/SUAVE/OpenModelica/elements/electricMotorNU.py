#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import ElectricMotorNUData
from SUAVE.modelica.src.hensup.elements.electricMotor import ElectricMotor
from SUAVE.modelica.src.hensup.elements.powerElectronics import PowerElectronics
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class ElectricMotorNU(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
        
        self.electricMotor = ElectricMotor('Electric_Motor_PU')
        self.powerElectronics = PowerElectronics('Power_Electronics_PU')
        
        self.unitSet = (self.electricMotor, self.powerElectronics)
        
        self.data = ElectricMotorNUData()
        
        self.buildParameters()
        
    def collectComponents(self, component):
        if component.id == 'parameter':
            self.parameterSet.add(component)
            self.parameterMoStringSet.add(component.moString)
        elif component.id == 'variable':
            self.variableSet.add(component)
            self.variableMoStringSet.add(component.moString)
        elif component.id == 'equation':
            self.equationMoStringSet.add(component.moString)
        elif component.id == 'init_equation':
            self.equationMoStringSet.add(component.moString)
        else:
            TypeError(f'Component type {component.typ} not supported.')
        
    def add_component(self, component):
        setattr(self, component.name, component)
        self.collectComponents(component)

        
    def setModelParameters(self, path=None):
        pass

    def buildVariables(self):
        pass
                
    def buildParameters(self):        
        self.add_component(Parameter('electricMotorGravimetricPower', 
                                     value=self.data.electricMotorGravimetricPower))
        self.add_component(Parameter('electricMotorEfficiency', 
                                     value=self.data.electricMotorEfficiency))
        
        self.add_component(Parameter('powerElectronicsGravimetricPower', 
                                     value=self.data.powerElectronicsGravimetricPower))
        self.add_component(Parameter('powerElectronicsEfficiency', 
                                     value=self.data.powerElectronicsEfficiency))
        
    def buildEquations(self):
        pass