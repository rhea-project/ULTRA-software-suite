# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:08:22 2020

@author: senth
"""

import math

def Lift_induced_drag(Cl , AR , e):
    
    Cdi = (Cl**2)/(math.pi*AR*e)
    
    return Cdi