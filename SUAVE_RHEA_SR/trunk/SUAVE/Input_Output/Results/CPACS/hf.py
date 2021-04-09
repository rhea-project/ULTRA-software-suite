# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 12:54:21 2019

@author: TUBS, ZJU
"""
import numpy as np

class HF():
    
  def Profiles(uID,Xcoordinates,Zcoordinates):
    XCoorStr = [str(i) for i in Xcoordinates]
    YCoorStr = [str(i) for i in np.zeros(np.shape(XCoorStr))]
    ZCoorStr = [str(i) for i in Zcoordinates]
    
    sep = ';'
    Str = "<wingAirfoil uID=\""+uID+"\">\n"\
"<name>NACA0009 Airfoil</name>\n"\
"<description>NACA 4-digit symmetrical airfoil with a thickness of 9percent</description>\n"\
"<pointList>\n"\
"<x mapType="+'"vector"'+">"+sep.join(XCoorStr)+"</x>\n"\
"<y mapType="+'"vector"'+">"+sep.join(YCoorStr)+"</y>\n"\
"<z mapType="+'"vector"'+">"+sep.join(ZCoorStr)+"</z>\n"\
"</pointList>\n"\
"</wingAirfoil>\n"

    return Str
  
  def Header():
    Str="<cpacs>\n"\
"<header>\n"\
"<name>genericSystem example file</name>\n"\
"<description>Example file for use of genericSystem components</description>\n"\
"<creator>Jonas Jepsen from TU-Hamburg-Harburg</creator>\n"\
"<timestamp>2017-05-31T08:55:23</timestamp>\n"\
"<version>1.0</version>\n"\
"<cpacsVersion>3.0</cpacsVersion>\n"\
"</header>\n"\
"<vehicles>\n"\
"<aircraft>\n"\
"<model uID="+'"Aircraft"'+">\n"\
"<name>genericSystem example</name>\n"\
"<description>Container to show use of genericSystem nodes</description>\n"
    
    return Str
    

       
  def Footer():    
       
    Str ="</wingAirfoils>\n"\
"<fuselageProfiles>\n"\
"<fuselageProfile uID="+'"fuselageCircle"'+">\n"\
"<name>Circle</name>\n"\
"<description>Profile build up from set of Points on Circle where may Dimensions are 1..-1</description>\n"\
"<pointList>\n"\
"<x mapType="+'"vector"'+">0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;</x>\n"\
"<y mapType="+'"vector"'+">0.5;0.49384;0.47553;0.4455;0.40451;0.35355;0.29389;0.227;0.15451;0.07822;0;-0.07822;-0.15451;-0.227;-0.29389;-0.35355;-0.40451;-0.4455;-0.47553;-0.49384;-0.5;-0.49384;-0.47553;-0.4455;-0.40451;-0.35355;-0.29389;-0.227;-0.15451;-0.07822;-0;0.07822;0.15451;0.227;0.29389;0.35355;0.40451;0.4455;0.47553;0.49384;0.5;</y>\n"\
"<z mapType="+'"vector"'+">0;0.07822;0.15451;0.227;0.29389;0.35355;0.40451;0.4455;0.47553;0.49384;0.5;0.49384;0.47553;0.4455;0.40451;0.35355;0.29389;0.227;0.15451;0.07822;0;-0.07822;-0.15451;-0.227;-0.29389;-0.35355;-0.40451;-0.4455;-0.47553;-0.49384;-0.5;-0.49384;-0.47553;-0.4455;-0.40451;-0.35355;-0.29389;-0.227;-0.15451;-0.07822;-0;</z>\n"\
"</pointList>\n"\
"</fuselageProfile>\n"\
"</fuselageProfiles>\n"\
"</profiles>\n"\
"</vehicles>\n"\
"</cpacs>\n"
    
    return Str
  
  def wingFooter():
     
   Str ="</model>\n"\
"</aircraft>\n"\
"<profiles>\n"\
"<wingAirfoils>\n"
     
   return Str