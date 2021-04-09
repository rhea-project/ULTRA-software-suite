#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import BatteryData
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class Battery(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
        
        self.data = BatteryData()
        
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
        
    def eqPowerOutput(self):
        return f'0 = {self.powerOutput.name} - ( {self.powerDischarge.name} - {self.powerCharge.name} )'
    
    def eqPowerDischarge(self):
        return f'0 = powerDischarge - powerStar * ( if powerStar < 0 then 1 else 0 )'
        
    def eqPowerCharge(self):
        return f'0 = powerCharge - powerStar * ( if powerStar >= 0 then -1 else 0 )'

    def buildVariables(self):
        self.add_component(Variable('powerStar'))
        self.add_component(Variable('powerOutput'))
        self.add_component(Variable('powerDischarge'))
        self.add_component(Variable('powerCharge'))
                
    def buildParameters(self):
        self.add_component(Parameter('capacityRated', 
                                     value=self.data.capacityRated))
        
        self.add_component(Parameter('gravimetricCapacity', 
                                     value=self.data.gravimetricCapacity))
        
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
        self.add_component(Equation('eqPowerOutput', self.eqPowerOutput()))
        self.add_component(Equation('eqPowerDischarge', self.eqPowerDischarge()))
        self.add_component(Equation('eqPowerCharge', self.eqPowerCharge()))