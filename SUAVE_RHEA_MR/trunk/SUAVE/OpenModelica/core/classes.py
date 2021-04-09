#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports

# related third party imports
    
# local application/library specific imports

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class Parameter(object):
    def __init__(self, name, typ='Real', comment='', value=None):
        self.id = 'parameter'
        self.name = name
        self.typ = typ
        self.comment = comment
        self.value = value
        
        self.moString = self.createMoString()
        
    def createMoString(self):
        return f'parameter {self.typ} {self.name} \"{self.comment}\";'

class Variable(object):
    def __init__(self, name, typ='Real', comment=''):
        self.id = 'variable'
        self.name = name
        self.typ = typ
        self.comment = comment  
        
        self.moString = self.createMoString()
        
    def createMoString(self):
        return f'{self.typ} {self.name} \"{self.comment}\";'
        
class Equation(object):
    def __init__(self, name, expr, comment=''):
        self.id = 'equation'
        self.name = name
        self.expr = expr
        self.comment = comment
        
        self.moString = self.createMoString()
        
    def createMoString(self):
        return f'{self.expr};'
    
class InitialEquation(Equation):
    def __init__(self, name, expr, comment=''):
        Equation.__init__(self, name, expr, comment=comment)
        self.id = 'init_equation'