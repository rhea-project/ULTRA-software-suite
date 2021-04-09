#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import DuctedFanData
from SUAVE.modelica.src.hensup.elements.expansionNozzle import ExpansionNozzle
from SUAVE.modelica.src.hensup.elements.compressionNozzle import CompressionNozzle
from SUAVE.modelica.src.hensup.elements.fan import Fan
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class DuctedFan(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
        
        self.inletNozzle = ExpansionNozzle('InletNozzle')
        self.fanNozzle = CompressionNozzle('FanNozzle')
        self.fan = Fan('Fan')
        self.data = DuctedFanData()
        
        self.unitSet = (self.inletNozzle, self.fanNozzle, self.fan)
        
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
        self.add_component(Parameter('thrustDesign', 
                                     value=self.data.thrustDesign))        
        self.add_component(Parameter('MachNumberDesign', 
                                     value=self.data.MachNumberDesign))
        self.add_component(Parameter('altitudeDesign', 
                                     value=self.data.altitudeDesign))
        
        self.add_component(Parameter('inletNozzleEfficiencyPoly', 
                                     value=self.inletNozzle.efficiencyPoly.value))        
        self.add_component(Parameter('inletNozzlePressureRatio', 
                                     value=self.inletNozzle.pressureRatio.value))
        
        self.add_component(Parameter('fanEfficiencyPoly', 
                                     value=self.fan.efficiencyPoly.value))
        self.add_component(Parameter('fanPressureRatio', 
                                     value=self.fan.pressureRatio.value))       
        
        self.add_component(Parameter('fanNozzleEfficiencyPoly', 
                                     value=self.fanNozzle.efficiencyPoly.value))
        self.add_component(Parameter('fanNozzleEfficiencyPoly', 
                                     value=self.fanNozzle.pressureRatio.value))    
    def buildEquations(self):
        pass