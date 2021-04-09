# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:48:53 2019

@author: TUBS, ZJU
"""

import numpy as np
from Attr import *

class Fuselage(attributes):
      
 def __init__(self,*args):   
    Str=""
    for fuse in args:
     self.uID=fuse[0]   
     self.Scaling=fuse[1]
     self.Rotation=fuse[2]
     self.Translation=fuse[3]
     self.Sections=fuse[4]
     self.Positioning=fuse[5]
     self.Segments=fuse[6]
     
     super().__init__(self.Scaling,self.Rotation,self.Translation)
     Str=Str+self.ConstructFuselage()
    Str = "<fuselages>\n"+Str+"</fuselages>\n"
    self.Str=Str
  
 def element(self,uID,airfoilUID,Sca,Rot,Tra):

    Str = "<element uID=\""+uID+"\">\n"\
"<name>"+uID+"</name>\n"\
"<description>"+uID+"</description>\n"\
"<profileUID>"+airfoilUID+"</profileUID>\n"+\
self.transformation(uID+"_Transf",Sca,Rot,Tra)+\
"</element>\n"    
    return Str

 def ConstructFuselage(self):
      
    Str = "<fuselage uID=\""+self.uID+"\">\n"\
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
    "</segments>\n"+\
"</fuselage>\n"
    return Str       
 
 def __str__(self):
     return self.Str