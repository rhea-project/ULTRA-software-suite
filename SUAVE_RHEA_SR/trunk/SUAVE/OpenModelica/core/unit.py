#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports
from hensup.src.core.template import Template

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class Unit(object):
    def __init__(self, name):
        self.name = name
        self.unitSet = set()
        self.parameterSet = set()
        self.variableSet = set()
        self.initEquationSet = set()
        self.equationSet = set()
        
        self.moFile = Template(self)
        
    def collectComponents(self, component):
        if component.id == 'parameter':
            self.parameterSet.add(component.moString)
        elif component.id == 'variable':
            self.variableSet.add(component.moString)
        elif component.id == 'equation':
            self.equationSet.add(component.moString)
        elif component.id == 'init_equation':
            self.equationSet.add(component.moString)
        else:
            TypeError(f'Component type {component.typ} not supported.')
            
    def add_component(self, component):
        setattr(self, component.name, component)
        self.collectComponents(component)

    def buildParameters(self):
        pass

    def buildVariables(self):
        pass

    def buildEquations(self):
        pass