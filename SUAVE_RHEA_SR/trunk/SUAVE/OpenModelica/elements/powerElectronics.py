#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import PowerElectronicsData
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class PowerElectronics(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
        
        self.data = PowerElectronicsData()
        
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

    def buildVariables(self):
        pass
                
    def buildParameters(self):
        self.add_component(Parameter('powerRated', 
                                     value=self.data.powerRated))
        
        self.add_component(Parameter('gravimetricPower', 
                                     value=self.data.gravimetricPower))
        self.add_component(Parameter('efficiency', 
                                     value=self.data.efficiency))
        
    def buildEquations(self):
        pass