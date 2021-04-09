#-----------------------------------------------------------------------------#
# Import modules
#-----------------------------------------------------------------------------#
# standard library imports
import os
import csv
import math
import collections

# related third party imports
import numpy as np

# local application/library specific imports
from src.hensup.database.constants import air, naturalConstants

#-----------------------------------------------------------------------------#
# Superclass for atmospheric models
#-----------------------------------------------------------------------------#
class AtmosModelMethods():
    def __init__(self):
        pass
    
    def calcTemperatureAtmos(self, T0, dTdh, h):
        """Calculates temperature in Troposphere."""
        return T0 - dTdh * h
        
    def calcPressureAtmos(self, p0, T0, dTdh, h):
        """Calculates pressure in Troposphere."""
        if dTdh == 0:
            return self.calcPressureStrato(p0, T0, h)
        else:
            return self.calcPressureTropo(p0, T0, dTdh, h)
    
    def calcDensityAtmos(self, rho0, T0, dTdh, h):
        """Calculates density in Troposphere."""
        if dTdh == 0:
            return self.calcDensityStrato(rho0, T0, h)
        else:
            return self.calcDensityTropo(rho0, T0, dTdh, h)
    
    def calcPressureTropo(self, p0, T0, dTdh, dh):
        """Calculates pressure in Stratosphere."""
        R = naturalConstants.gasConstant / air.molarMass
        return p0 * (1 - dTdh * dh / T0)**(naturalConstants.standardAcceleration / (R * dTdh))
    
    def calcDensityTropo(self, rho0, T0, dTdh, dh):
        """Calculates density in Stratosphere."""
        R = naturalConstants.gasConstant / air.molarMass
        return rho0 * (1 - dTdh * dh / T0)**(naturalConstants.standardAcceleration / (R * dTdh) - 1)
    
    def calcPressureStrato(self, p0, T0, dh):
        """Calculates pressure in Stratosphere."""
        R = naturalConstants.gasConstant / air.molarMass
        return p0 * math.exp(-naturalConstants.standardAcceleration / (R * T0) * dh)
    
    def calcDensityStrato(self, rho0, T0, dh):
        """Calculates density in Stratosphere."""
        R = naturalConstants.gasConstant / air.molarMass
        return rho0 * math.exp(-naturalConstants.standardAcceleration / (R * T0) * dh)  
        
#-----------------------------------------------------------------------------#
# International Standard Atmosphere (ISA)
#-----------------------------------------------------------------------------#
class ISA(AtmosModelMethods):
    def __init__(self):
        self.groundTemperature = 288.15
        self.temperatureGradient = 0.0065
        self.heightTropo = 11000
            
    def calcISATemperature(self, height):   
        """Returns temperature of ISA model."""       
        ambientTemperature = []
        
        if not isinstance(height, collections.Iterable):
            height = [height]
            
        for h in height:
            T0 = self.groundTemperature
            dTdh = self.temperatureGradient
            
            if h <= self.heightTropo:
                dh = h
            else:
                dh = self.heightTropo
            
            T = self.calcTemperatureAtmos(T0, dTdh, dh)
            
            ambientTemperature.append(T)
    
        return np.asarray(ambientTemperature)
            
    
    def calcISAPressure(self, height):    
        """Returns pressure of ISA model."""      
        ambientPressure = []
        
        if not isinstance(height, collections.Iterable):
            height = [height]
            
        for h in height:
            p0 = self.groundPressure
            dTdh = self.temperatureGradient
            
            if h <= self.heightTropo:
                T0 = self.groundTemperature
                dh = h
            else:
                T0 = self.calcTemperature(self.heightTropo)
                dh = self.heightTropo
            
            p = self.calcPressureAtmos(p0, T0, dTdh, dh)
            
            ambientPressure.append(p)
    
        return np.asarray(ambientPressure)
    
    
    def calcISADensity(self, height):   
        """Returns density of ISA model."""       
        ambientDensity = []
        
        if not isinstance(height, collections.Iterable):
            height = [height]
            
        for h in height:
            rho0 = self.groundDensity
            dTdh = self.temperatureGradient
            
            if h <= self.heightTropo:
                T0 = self.groundTemperature
                dh = h
            else:
                T0 = self.calcTemperature(self.heightTropo)
                dh = self.heightTropo
            
            rho = self.calcDensityAtmos(rho0, T0, dTdh, dh)
                
            ambientDensity.append(rho)
    
        return np.asarray(ambientDensity)
    
   
#-----------------------------------------------------------------------------#
# Artic Minimum Standard Atmosphere (AMSA)
#-----------------------------------------------------------------------------#    
class AMSA(AtmosModelMethods):
    def __init__(self):
        self.groundTemperature = 223.15
        self.temperatureGradient = [-0.01, 0, 0.00472]
        self.heightTropo = [1500, 3000, 15500]
    
    def calcAMSATemperature(self, height): 
        """Returns temperature of AMSA model."""       
        ambientTemperature = []     
        
        
        if not isinstance(height, collections.Iterable):
            height = [height]
            
        for h in height:
            if h <= self.heightTropo[0]:
                T0 = self.groundTemperature
                dTdh = self.temperatureGradient[0]
                dh = h
                
            elif h > self.heightTropo[0] and h <= self.heightTropo[1]:
                T0 = self.calcTemperature(self.heightTropo[0])
                dTdh = self.temperatureGradient[1]
                dh = h - self.heightTropo[0]
                
            elif h > self.heightTropo[1] and h <= self.heightTropo[2]:
                T0 = self.calcTemperature(self.heightTropo[1])
                dTdh = self.temperatureGradient[2]
                dh = h - self.heightTropo[1]
                
            else: 
                T0 = self.calcTemperature(self.heightTropo[2])
                dTdh = self.temperatureGradient[2]
                dh = self.heightTropo[2] - self.heightTropo[1]
                
            T = self.calcTemperatureAtmos(T0, dTdh, dh)
            ambientTemperature.append(T)
            
        return np.asarray(ambientTemperature)

        
    def calcAMSAPressure(self, height):      
        """Returns pressure of AMSA model."""  
        ambientPressure = []     
        
        if not isinstance(height, collections.Iterable):
            height = [height]
                        
        for h in height:
            if h <= self.heightTropo[0]:
                p0 = self.groundPressure
                T0 = self.groundTemperature
                dTdh = self.temperatureGradient[0]
                dh = h
        
            elif h > self.heightTropo[0] and h <= self.heightTropo[1]:
                p0 = self.calcPressure(self.heightTropo[0])
                T0 = self.calcTemperature(self.heightTropo[0])
                dh = h - self.heightTropo[0]
                
            elif h > self.heightTropo[1] and h <= self.heightTropo[2]:
                p0 = self.calcPressure(self.heightTropo[1])
                T0 = self.calcTemperature(self.heightTropo[1])
                dTdh = self.temperatureGradient[2]
                dh = h - self.heightTropo[1]
                
            else:
                p0 = self.calcPressure(self.heightTropo[2])
                T0 = self.calcTemperature(self.heightTropo[2])
                dh = self.heightTropo[2] - self.heightTropo[1]
                
            p = self.calcPressureAtmos(p0, T0, dTdh, dh)
            
            ambientPressure.append(p)
    
        return np.asarray(ambientPressure)
                    
    
    def calcAMSADensity(self, height):   
        """Returns density of AMSA model."""     
        ambientDensity = []
        
        if not isinstance(height, collections.Iterable):
            height = [height]
            
        for h in height:
            if h <= self.heightTropo[0]:
                rho0 = self.groundDensity
                T0 = self.groundTemperature
                dTdh = self.temperatureGradient[0]
                dh = h
        
            elif h > self.heightTropo[0] and h <= self.heightTropo[1]:
                rho0 = self.calcDensity(self.heightTropo[0])
                T0 = self.calcTemperature(self.heightTropo[0])
                dh = h - self.heightTropo[0]
                
            elif h > self.heightTropo[1] and h <= self.heightTropo[2]:
                rho0 = self.calcDensity(self.heightTropo[1])
                T0 = self.calcTemperature(self.heightTropo[1])
                dTdh = self.temperatureGradient[2]
                dh = h - self.heightTropo[1]
                
            else:
                rho0 = self.calcDensity(self.heightTropo[2])
                T0 = self.calcTemperature(self.heightTropo[2])
                dh = self.heightTropo[2] - self.heightTropo[1]
            
            rho = self.calcDensityAtmos(rho0, T0, dTdh, dh)    
            ambientDensity.append(rho)
    
        return np.asarray(ambientDensity)
        
            
#-----------------------------------------------------------------------------#
# Tropical Maximum Standard Atmosphere (TMSA)
#-----------------------------------------------------------------------------#
class TMSA(AtmosModelMethods):
    def __init__(self):
        self.groundTemperature = 318.15
        self.temperatureGradient = 0.0065
        self.heightTropo = 11540        
        
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class AtmosModel(ISA, AMSA, TMSA):
    def __init__(self, atmosModel):
        self.groundPressure = 101325
        self.groundDensity = 1.225 
        
        if atmosModel == 0:
            AMSA.__init__(self)
            self.funcCalcTemperature = self.calcAMSATemperature
            self.funcCalcPressure = self.calcAMSAPressure
            self.funcCalcDensity = self.calcAMSADensity
        elif atmosModel == 1:
            ISA.__init__(self)
            self.funcCalcTemperature = self.calcISATemperature
            self.funcCalcPressure = self.calcISAPressure
            self.funcCalcDensity = self.calcISADensity
        elif atmosModel == 2:
            TMSA.__init__(self)
            self.funcCalcTemperature = self.calcISATemperature
            self.funcCalcPressure = self.calcISAPressure
            self.funcCalcDensity = self.calcISADensity
            
        self.atmosModel = atmosModel
        
        
    def calcTemperature(self, height):
        """Returns temperature dependent on the selected atmospheric model."""
        return self.funcCalcTemperature(height)

    def calcPressure(self, height):
        """Returns pressure dependent on the selected atmospheric model."""
        return self.funcCalcPressure(height)

    def calcDensity(self, height):
        """Returns density dependent on the selected atmospheric model."""
        return self.funcCalcDensity(height)


    def buildAtmosphericProperties(self):
        """Creates mission profile.
        
        buildAtmosphericProperties(self) ambient temperature, pressure and 
        density to the object.
        
        Parameters
        ----------
        self : `object`
            Instance of class.
        
        Returns
        -------
        self : `object`
            Object passed as input with additional attributes:
          
                - ``.ambientTemperature``: contains ambient temperature at 
                mission time steps (`numpy.ndarray`).
          
                - ``.ambientPressure``: contains ambient pressure at 
                mission time steps (`numpy.ndarray`).
          
                - ``.ambientDensity``: contains ambient density at 
                mission time steps (`numpy.ndarray`).
         
        Notes
        -----        
        
        """
        
        HStrato = 16000
    
        if max(self.height) > HStrato:
            raise ValueError('No atmospheric model available for altitudes \
                             greater than {}!'.format(HStrato))
                    
        self.ambientDensity = self.calcDensity(self.height)
        self.ambientPressure = self.calcPressure(self.height) / 1e+05
        self.ambientTemperature = self.calcTemperature(self.height)
        
