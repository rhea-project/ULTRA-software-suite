# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:03:20 2019

@author: TUBS, ZJU
"""

class attributes():
    
 def __init__(self,Scaling,Rotation,Translation):
     self.Scaling=Scaling
     self.Rotation=Rotation
     self.Traslation=Translation
 
 def transformation(self,uID,Scaling,Rotation,Translation):
    
    Str = "<transformation uID=\""+uID+"\">\n"\
"<scaling uID=\""+uID+"_Sca\">\n"\
"<x>"+str(Scaling[0])+"</x>\n"\
"<y>"+str(Scaling[1])+"</y>\n"\
"<z>"+str(Scaling[2])+"</z>\n"\
"</scaling>\n"\
"<rotation uID=\""+uID+"_Rot\">\n"\
"<x>"+str(Rotation[0])+"</x>\n"\
"<y>"+str(Rotation[1])+"</y>\n"\
"<z>"+str(Rotation[2])+"</z>\n"\
"</rotation>\n"\
"<translation uID=\""+uID+"_Tra\">\n"\
"<x>"+str(Translation[0])+"</x>\n"\
"<y>"+str(Translation[1])+"</y>\n"\
"<z>"+str(Translation[2])+"</z>\n"\
"</translation>\n"\
"</transformation>\n"
    
    
    return Str


 def element(self,uID,airfoilUID,Scaling,Rotation,Translation):
    
    Str = "<element uID=\""+uID+"\">\n"\
"<name>"+uID+"</name>\n"\
"<description>"+uID+"</description>\n"\
"<airfoilUID>"+airfoilUID+"</airfoilUID>\n"+\
self.transformation(uID+"_Transf",Scaling,Rotation,Translation)+\
"</element>\n"
     
    return Str

 def positioning(self,uID,length,sweep,dihedral,fromSection,toSection):
    
    Str = "<positioning uID=\""+uID+"\">\n"\
"<name>"+uID+"</name>\n"\
"<description>"+uID+"</description>\n"\
"<length>"+str(length)+"</length>\n"\
"<sweepAngle>"+str(sweep)+"</sweepAngle>\n"\
"<dihedralAngle>"+str(dihedral)+"</dihedralAngle>\n"\
"<fromSectionUID>"+fromSection+"</fromSectionUID>\n"\
"<toSectionUID>"+toSection+"</toSectionUID>\n"\
"</positioning>\n"
     
    return Str

 def segment(self,uID,fromSection,toSection):
    
    Str = "<segment uID=\""+uID+"\">\n"\
"<name>"+uID+"</name>\n"\
"<description>"+uID+"</description>\n"\
"<fromElementUID>"+fromSection+"</fromElementUID>\n"\
"<toElementUID>"+toSection+"</toElementUID>\n"\
"</segment>\n"
    
    return Str

 def section(self,uID,airfoiluID,Scaling,Rotation,Translation,elementScaling,elementRotation,elementTranslation):
    
    Str = "<section uID=\""+uID+"\">\n"\
"<name>"+uID+"</name>\n"\
"<description>"+uID+"</description>\n"+\
self.transformation(uID+"_Transf",Scaling,Rotation,Translation)+\
"<elements>\n"+\
self.element(uID+"_element",airfoiluID,elementScaling,elementRotation,elementTranslation)+\
"</elements>\n"\
"</section>\n"
    
    return Str