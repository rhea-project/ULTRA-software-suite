#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class MoTemplate(object):
    def __init__(self, unit):
        self.unit = unit
    
    def buildMoText(self):
        unit = self.unit
        text = []
        start = f'model {unit.name} \n'
        text.append(start)
        
        text.append(f'\n')
        text.append(f'\t\\\ Load sub-units \n')
        for unit in unit.unitSet:
            text.append(f'\t{unit.name} {unit.name}\n')
            
        text.append(f'\n')
        text.append(f'\t\\\ Declaring model parameter \n')
        for param in unit.parameterSet:
            text.append(f'\t{param} \n')
        
        text.append(f'\n')
        text.append(f'\t\\\ Declaring model variables \n')
        for var in unit.variableSet:
            text.append(f'\t{var} \n')
        
        text.append(f'\n')
        text.append(f'\t\\\ Declaring initial model equations \n')
        text.append(f'\tinitial equation \n')
        for eq in unit.initEquationSet:
            text.append(f'\t\t{eq} \n')
            
        text.append(f'\n')
        text.append(f'\t\\\ Declaring model equations \n')
        text.append(f'\tequation \n')
        for eq in unit.equationSet:
            text.append(f'\t\t{eq} \n')
        
        text.append(f'\n')
        end = f'end {unit.name};'
        text.append(end)
        
        self.text = text