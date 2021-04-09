# -*- coding: utf-8 -*-
"""
From SUAVE to CPACS


@author: TUBS, ZJU
"""

from Wing import *
import pandas as pd
from suave_demo import *
from hf import *
from Fuselage import *
import numpy as np

Aircraft= vehicle_setup()

WingProfile = pd.read_csv('HL.txt',header = 0,delim_whitespace=True)
VSProfile = pd.read_csv('NACA.txt',header = 0,delim_whitespace=True)
AircraftTag = 'Aircraft'

'''
Function to extract parameters
'''
def ArrangeParameters(Segmentation,Span,Root,Tip,Sweep,Dihedral,Wingtag):

 if (Segmentation):
  Sections = []
  Positions = [["position_1",0.0,0.0,0.0,"",WingTag+"_section_1"]]
  Segments = []
  for i,Val in enumerate(Segmentation.values()):
    Scaling = Root*Val.root_chord_percent
    Element = [[Scaling,Scaling,Scaling],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    temp1 = [Val.tag,"HighLift",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    temp1.extend(Element)
    Sections.append(temp1)
    if (i>0):
     SectionA =  WingTag+"_section_"+str(i)
     SectionB =WingTag+"_section_"+str(i+1)
     Length = Span*((Segmentation.values()[i].percent_span_location-Segmentation.values()[i-1].percent_span_location))/2
     Sweep = np.rad2deg(Val.sweeps.quarter_chord)
     Dihedral = np.rad2deg(Val.dihedral_outboard)
     temp2 = ["position_"+str(i+1)]
     temp2.extend([Length,Sweep,Dihedral,SectionA,SectionB])
     Positions.append(temp2)
     Segments.append(["segment_"+str(i),SectionA+"_element",SectionB+"_element"])

 else:
  RootScaling = Root   
  TipScaling = Tip
  Length = Span/2
  Element1 = [[RootScaling,RootScaling,RootScaling],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  Element2 = [[TipScaling,TipScaling,TipScaling],[0.0,0.0,0.0],[0.0,0.0,0.0]]
 
  Section1 = "section_1"
  Section2 = "section_2"
 
  temp1 = [Section1,"HighLift",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]   
  temp2 = [Section2,"HighLift",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
  temp1.extend(Element1)
  temp2.extend(Element2)
  Sections = [temp1,temp2]

  PosA = ["position_1",0.0,0.0,0.0,""]
  PosA.extend([WingTag+"_"+Section1])
  PosB = ["position_2",Length,np.rad2deg(Sweep),Dihedral]
  PosB.extend([WingTag+"_"+Section1,WingTag+"_"+Section2])
  Positions = [PosA,PosB]
  Segments =[["segment_1",WingTag+"_"+Section1+"_element",WingTag+"_"+Section2+"_element"]]
 
 return Sections, Positions, Segments

'''
Function to convert quarter chord to leading edge sweep
'''

def qcleSweep(sweep,span,root,tip):
    
    chord_fraction      = 0.25
    newsweep = np.arctan(((root*chord_fraction)+(np.tan(sweep)*span)-(tip*chord_fraction))/span)
    
    return newsweep
    
    

'''''
Wing Details
1. Sectioning
2. Elements
3. Positioning
4. Segments
'''''

#Main Wing
Origin = Aircraft.wings.main_wing.origin
Segmentation = Aircraft.wings.main_wing.Segments
Root = Aircraft.wings.main_wing.chords.root
Tip = Aircraft.wings.main_wing.chords.tip
Span = Aircraft.wings.main_wing.spans.projected
LESweep = qcleSweep(Aircraft.wings.main_wing.sweeps.quarter_chord,Span,Root,Tip)
Sweep = Aircraft.wings.main_wing.sweeps.quarter_chord
Span = Aircraft.wings.main_wing.spans.projected/np.cos(LESweep)
Dihedral = 5
WingTag =  AircraftTag+'_'+Aircraft.wings.main_wing.tag
Symmetry = True
Sections, Positions, Segments = ArrangeParameters(Segmentation,Span,Root,Tip,LESweep,Dihedral,WingTag)
MainWing = [WingTag,[1.0,1.0,1.0],[0.0,0.0,0.0],Origin,Sections, Positions, Segments,Symmetry]

#Horizontal Stabilizer
Origin = Aircraft.wings.horizontal_stabilizer.origin
Segmentation = Aircraft.wings.horizontal_stabilizer.Segments
Root = Aircraft.wings.horizontal_stabilizer.chords.root
Tip = Aircraft.wings.horizontal_stabilizer.chords.tip
Sweep = Aircraft.wings.horizontal_stabilizer.sweeps.quarter_chord
Span = Aircraft.wings.horizontal_stabilizer.spans.projected
Dihedral = 5
WingTag =  AircraftTag+'_'+Aircraft.wings.horizontal_stabilizer.tag
Symmetry = True
Sections, Positions, Segments = ArrangeParameters(Segmentation,Span,Root,Tip,Sweep,Dihedral,WingTag)
HS = [WingTag,[1.0,1.0,1.0],[0.0,0.0,0.0],Origin,Sections, Positions, Segments,Symmetry]

#Vertical Stabilizer
Origin = Aircraft.wings.vertical_stabilizer.origin
Segmentation = Aircraft.wings.vertical_stabilizer.Segments
Root = Aircraft.wings.vertical_stabilizer.chords.root
Tip = Aircraft.wings.vertical_stabilizer.chords.tip
Sweep = Aircraft.wings.vertical_stabilizer.sweeps.quarter_chord
Span = Aircraft.wings.vertical_stabilizer.spans.projected
Dihedral = 0
WingTag =  AircraftTag+'_'+Aircraft.wings.vertical_stabilizer.tag
Symmetry = False
Sections, Positions, Segments = ArrangeParameters(Segmentation,Span,Root,Tip,Sweep,Dihedral,WingTag)
VS = [WingTag,[1.0,1.0,1.0],[90.0,0.0,0.0],Origin,Sections, Positions, Segments,Symmetry]

#Genrate xml for wing
WingXml = Wing(MainWing,HS,VS).__str__()


FuselageTag =  AircraftTag+'_'+Aircraft.fuselages.fuselage.tag
NoseSec = []
NosePos = [["position_1",0.0,90.0,0.0,"",FuselageTag+"_section_1"]]
NoseSeg =[]

for i,Val in enumerate(NoseScaling):
    NoseSec.append(["section_"+str(i+1),"fuselageCircle",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[Val,Val,Val],[0.0,0.0,0.0],[0.0,0.0,NoseCenterz[i]]])
for i,Val in enumerate(NoseCenterx,start = 1):
    NosePos.append(["position_"+str(i+1),Val,90,0.0,FuselageTag+"_section_"+str(i),FuselageTag+"_section_"+str(i+1)])
for i, Val in enumerate(NoseCenterx,start=1):
    NoseSeg.append(["segment_"+str(i),FuselageTag+"_section_"+str(i)+"_element",FuselageTag+"_section_"+str(i+1)+"_element"])

for i,Val in enumerate(CabinScaling,start = N+1):
    NoseSec.append(["section_"+str(i),"fuselageCircle",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[Val,Val,Val],[0.0,0.0,0.0],[0.0,0.0,0.0]])
for i,Val in enumerate(CabinCenterx,start = N):
    NosePos.append(["position_"+str(i),Val,90,0.0,FuselageTag+"_section_"+str(i-1),FuselageTag+"_section_"+str(i)])
for i,Val in enumerate(CabinCenterx,start = N-1):
    NoseSeg.append(["segment_"+str(i),FuselageTag+"_section_"+str(i)+"_element",FuselageTag+"_section_"+str(i+1)+"_element"])
#    
for i,Val in enumerate(TailScaling):
    NoseSec.append(["section_"+str(i+N+2),"fuselageCircle",[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[Val,Val,Val],[0.0,0.0,0.0],[0.0,0.0,TailCenterz[i]]])
for i,Val in enumerate(TailCenterx,start = N+1):
    NosePos.append(["position_"+str(i),Val,90,0.0,FuselageTag+"_section_"+str(i-1),FuselageTag+"_section_"+str(i)])
for i,Val in enumerate(TailCenterx,start = N):
    NoseSeg.append(["segment_"+str(i),FuselageTag+"_section_"+str(i)+"_element",FuselageTag+"_section_"+str(i+1)+"_element"])

Fuse = [FuselageTag,[1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],NoseSec,NosePos,NoseSeg]

FuselageXml = Fuselage(Fuse).__str__()

xml_tmp = HF.Header() + WingXml+FuselageXml + HF.wingFooter() + HF.Profiles('HighLift',WingProfile['Xcoor'],WingProfile['Ycoor']) + HF.Profiles('NACA',VSProfile['Xcoor'],VSProfile['Ycoor'])+HF.Footer()


import xml.dom.minidom

uglyxml = xml_tmp

xml = xml.dom.minidom.parseString(uglyxml)

fpretty = xml.toprettyxml(newl='') # xml_pretty_str


f= open("suave_2_cpacs_demo.xml","w")

f.write(fpretty)
f.close()
