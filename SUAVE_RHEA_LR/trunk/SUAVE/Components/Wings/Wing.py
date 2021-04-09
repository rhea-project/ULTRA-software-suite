## @ingroup Components-Wings
# Wing.py
# 
# Created:  
# Modified: Sep 2016, E. Botero
#           Jul 2017, M. Clarke
#           Oct 2017, E. Botero
#           Oct 2018, T. MacDonald
#           Oct 2019, Y. Liu
#           Dec 2019, S. Karpuk

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import SUAVE
from SUAVE.Core import Data
from SUAVE.Components import Component, Lofted_Body, Mass_Properties
from .Airfoils import Airfoil

# ------------------------------------------------------------
#   Wing
# ------------------------------------------------------------

## @ingroup Components-Wings
class Wing(Lofted_Body):
    """This class defines the wing in SUAVE

    Assumptions:
    None

    Source:
    N/A

    Inputs:
    None

    Outputs:
    None

    Properties Used:
    N/A
    """      
    def __defaults__(self):
        """This sets the default values of a wing defined in SUAVE.
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """         

        self.tag             = 'wing'
        self.mass_properties = Mass_Properties()
        self.origin          = [[0.0,0.0,0.0]]

        self.number_of_wings = 1

        self.EMWET                     = False
        
        self.symmetric                 = True
        self.vertical                  = False
        self.t_tail                    = False
        
        self.wing_loading              = 0.0
        self.GLOV                      = 0.0            # Wing glove area (sq m)           
        self.FSTRT                     = 0.0            # Wing strut bracing factor (0.0 - nowing strut to
                                                        #                            1.0 - full benefit from strut bracing
        self.FAERT                     = 0.0            # Aeroelastic tailoring factor used in the design of the wing (0.0 -  no aeroelastic tailoring to
                                                        #                                                              1.0 -  maximum aeroelastic tailoring)
        self.PCTL                      = 0.5            # Percent of load carried by the wing
        self.composite_fraction        = 0.0
        self.taper                     = 0.0
        self.dihedral                  = 0.0
        self.aspect_ratio              = 0.0
        self.thickness_to_chord        = 0.0
        self.span_efficiency           = 0.9
        self.aerodynamic_center        = [0.0,0.0,0.0]
        self.exposed_root_chord_offset = 0.0

        self.spans = Data()
        self.spans.projected = 0.0
        
        # Define areas
        self.areas = Data()
        self.areas.reference = 0.0
        self.areas.exposed   = 0.0
        self.areas.affected  = 0.0
        self.areas.wetted    = 0.0
        self.areas.flap      = 0.0

        # Define chords
        self.chords = Data()
        self.chords.mean_aerodynamic = 0.0
        self.chords.mean_geometric   = 0.0
        self.chords.root             = 0.0
        self.chords.tip              = 0.0
        
        # Define sweeps
        self.sweeps               = Data()
        self.sweeps.quarter_chord = 0.0
        self.sweeps.leading_edge  = 0.0
        self.sweeps.half_chord    = 0.0
        self.sweeps.variable       = 0                      # Variable sweep indicator (1 for variable sweepwing
                                                            #                           0 for fixed wing)

        # Define twists
        self.twists = Data()
        self.twists.root = 0.0
        self.twists.tip  = 0.0

        # Define flaps
        self.control_surfaces = Data()
        self.flaps = Data()
        self.flaps.chord      = 0.0
        self.flaps.angle      = 0.0
        self.flaps.span_start = 0.0
        self.flaps.span_end   = 0.0
        self.flaps.type       = None
        self.flaps.area       = 0.0

        # Define slats
        self.slats = Data()
        self.slats.chord      = 0.0
        self.slats.angle      = 0.0
        self.slats.span_start = 0.0
        self.slats.span_end   = 0.0
        self.slats.type       = None

        self.high_lift     = False
        self.high_mach     = False
        self.vortex_lift   = False

        # Define wing folding characteristics    
        self.folding       = False
        self.folding_penalty = 0.0 

        # Define wing transition
        self.transition_x_upper = 0.0
        self.transition_x_lower = 0.0
        
        self.Airfoil            = Data()
        self.Segments           = SUAVE.Core.ContainerOrdered()
        self.Fuel_Tanks         = SUAVE.Core.Container()

        # Defined for the vertical tail
        self.wing_mount_coef    = 0                                 # Wing mount coefficient (0 - conventional tail to
                                                                    #                         1 - T-tail)
        
        # defined for BWB outerwing
        self.outerwing = Data()
        self.outerwing.spans = 0.0 
        self.outerwing.thickness_to_chord = 0.0
        self.outerwing.sweeps_quarter_chord = 0.0
        self.outerwing.taper = 0.0 
        self.outerwing.chords_root = 0.0
        self.outerwing.chords_tip = 0.0
        self.outerwing.chords_mean_aerodynamic = 0.0
        
        # defined for BWB aft center body        
        self.aft_body = Data()
        self.aft_centerbody_area               = 0.0
        self.aft_centerbody_taper              = 0.0
        self.aft_body.aft_spar_location_center = 0.7
        self.aft_body.aft_spar_location_side   = 0.7
        self.aft_body.spans                    = 0.0 
        self.aft_body.thickness_to_chord       = 0.0
        self.aft_body.sweeps_quarter_chord     = 0.0
        self.aft_body.taper                    = 0.0
        
        # defined for BWB cabin        
        self.cabin = Data()
        self.cabin.spans                = 0.0
        self.cabin.number_of_bays       = 0       # Number of bays for BWB
        self.cabin.thickness_to_chord   = 0.0
        self.cabin.sweeps_quarter_chord = 0.0
        self.cabin.taper                = 0.0
        self.cabin_sweep_leading_edge   = 0.0     # Cabin sweep angle (rad)
        self.cabin_length               = 0.0

        # Define wing material properties (for EMWET)
        self.materials = Data()
        self.materials.composites = False

        # Define aeroelastic inputs
        self.aeroelasticity = Data()
        self.aeroelasticity.CG_offset = 0.0


    def append_segment(self,segment):
        """ Adds a segment to the wing 
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """ 

        # Assert database type
        if not isinstance(segment,Data):
            raise Exception('input component must be of type Data()')

        # Store data
        self.Segments.append(segment)

        return
    
    def append_airfoil(self,airfoil):
        """ Adds an airfoil to the segment 
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """ 

        # Assert database type
        if not isinstance(airfoil,Data):
            raise Exception('input component must be of type Data()')

        # Store data
        self.Airfoil.append(airfoil)

        return        


    def append_control_surface(self,control_surface):
        """ Adds a component to vehicle 
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """ 

        # Assert database type
        if not isinstance(control_surface,Data):
            raise Exception('input control surface must be of type Data()')

        # Store data
        self.control_surfaces.append(control_surface)

        return
    
    def append_fuel_tank(self,fuel_tank):
        """ Adds a fuel tank to the wing 
    
        Assumptions:
        None

        Source:
        N/A

        Inputs:
        None

        Outputs:
        None

        Properties Used:
        N/A
        """ 

        # Assert database type
        if not isinstance(fuel_tank,Data):
            raise Exception('input component must be of type Data()')

        # Store data
        self.Fuel_Tanks.append(fuel_tank)

        return
