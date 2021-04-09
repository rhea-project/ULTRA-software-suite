# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:15:18 2019

@author: Uday Velakur
"""
import numpy as np
from Attr import *

class Wing(attributes):
       
 def __init__(self,*args):   
     Str = "" 
     for wing in args:
      self.uID=wing[0]   
      self.Scaling=wing[1]
      self.Rotation=wing[2]
      self.Translation=wing[3]
      self.Sections=wing[4]
      self.Positioning=wing[5]
      self.Segments=wing[6]
      self.Symmetry=wing[7]
     
      super().__init__(self.Scaling,self.Rotation,self.Translation)
      Str=Str+self.ConstructWing()
     Str = "<wings>\n"+Str+"</wings>\n"
     self.Str=Str
 
 def ConstructWing(self):

    Sym = "  symmetry="+'"x-z-plane"' if(self.Symmetry) else ""
    Str ="<wing uID=\""+self.uID+"\""+Sym+">\n"\
"<name>"+self.uID+"</name>\n"\
"<description>"+self.uID+"</description>\n"+\
self.transformation(self.uID+"_Transf",self.Scaling,self.Rotation,self.Traslation)+\
"<sections>\n"
 
    temp1 = ""
    temp2 = ""
    temp3 = ""
    
    for _,val in enumerate(self.Sections):
        temp1 = temp1+\
        self.section(self.uID+'_'+val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7])        
     
    Str = Str+temp1+\
    "</sections>\n"+\
    "<positionings>\n"
    
    for _, val in enumerate(self.Positioning):
        temp2 = temp2+\
        self.positioning(self.uID+'_'+val[0],val[1],val[2],val[3],val[4],val[5])

    Str = Str+temp2+\
    "</positionings>\n"+\
    "<segments>\n"
    
    for _, val in enumerate(self.Segments):
        temp3 = temp3+\
        self.segment(self.uID+'_'+val[0],val[1],val[2])

    Str = Str+temp3+\
"</segments>\n"+"</wing>\n"

    return Str
 
 def __str__(self):
     return self.Str