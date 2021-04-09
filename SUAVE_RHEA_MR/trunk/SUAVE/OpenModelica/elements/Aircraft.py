#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from SUAVE.modelica.data.component_data import AircraftData
from SUAVE.modelica.src.hensup.core.classes import Parameter, Variable, Equation

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class Aircraft(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.parameterMoStringSet = set()
        self.variableSet = set()
        self.variableMoStringSet = set()
        self.initEquationMoStringSet = set()
        self.equationMoStringSet = set()
            
        self.unitSet = ()
        
        self.data = AircraftData()
        
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
        self.add_component(Parameter('keroseneTankCapacityRated', 
                                     value=self.data.keroseneTankCapacityRated))
        self.add_component(Parameter('batteryBatteryNUCapacityRated', 
                                     value=self.data.batteryBatteryNUCapacityRated))
        
        self.add_component(Parameter('powerElectronicsBatteryNUPowerRated', 
                                     value=self.data.powerElectronicsBatteryNUPowerRated))
        self.add_component(Parameter('thrustTurboFanDesignPoint', 
                                     value=self.data.thrustTurboFanDesignPoint))
        self.add_component(Parameter('thrustDuctedFanDesignPoint', 
                                     value=self.data.thrustDuctedFanDesignPoint))
        
        self.add_component(Parameter('electricMotorElectricMotorNUPowerRated', 
                                     value=self.data.electricMotorElectricMotorNUPowerRated))
        self.add_component(Parameter('powerElectronicsElectricMotorNUPowerRated', 
                                     value=self.data.powerElectronicsElectricMotorNUPowerRated))
        
        self.add_component(Parameter('bypassRatio', 
                                     value=self.data.bypassRatio))
        self.add_component(Parameter('compressorPressureRatio', 
                                     value=self.data.compressorPressureRatio))
        self.add_component(Parameter('temperatureInletTurbine', 
                                     value=self.data.temperatureInletTurbine))
        
    def buildEquations(self):
        pass