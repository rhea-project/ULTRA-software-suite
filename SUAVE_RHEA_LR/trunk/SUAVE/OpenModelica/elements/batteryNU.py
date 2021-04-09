#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import BatteryNUData
from SUAVE.modelica.src.hensup.elements.battery import Battery
from SUAVE.modelica.src.hensup.elements.powerElectronics import PowerElectronics
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class BatteryNU(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
        
        self.battery = Battery('Battery_PU')
        self.powerElectronics = PowerElectronics('Power_Electronics_PU')
            
        self.unitSet = (self.battery, self.powerElectronics)
        
        self.data = BatteryNUData()
        
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
        self.add_component(Parameter('powerElectronicsGravimetricPower', 
                                     value=self.data.powerElectronicsGravimetricPower))
        self.add_component(Parameter('powerElectronicsEfficiency', 
                                     value=self.data.powerElectronicsEfficiency))
        
        self.add_component(Parameter('selfDischargeRate', 
                                     value=self.data.selfDischargeRate))
        self.add_component(Parameter('specificPowerCharge', 
                                     value=self.data.specificPowerCharge))
        self.add_component(Parameter('specificPowerDischarge', 
                                     value=self.data.specificPowerDischarge))
        
        self.add_component(Parameter('stateOfEnergyMax', 
                                     value=self.data.stateOfEnergyMax))
        self.add_component(Parameter('stateOfEnergyMin', 
                                     value=self.data.stateOfEnergyMin))
        self.add_component(Parameter('stateOfEnergyDelta', 
                                     value=self.data.stateOfEnergyDelta))
        self.add_component(Parameter('stateOfEnergyInit', 
                                     value=self.data.stateOfEnergyInit))
        
        self.add_component(Parameter('efficiencyCharge', 
                                     value=self.data.efficiencyCharge))
        self.add_component(Parameter('efficiencyDischarge', 
                                     value=self.data.efficiencyDischarge))
        
    def buildEquations(self):
        pass