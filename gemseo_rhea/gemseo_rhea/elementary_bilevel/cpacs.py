

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSStructureData


class CPACS_MR_SBW(CPACSStructureData):

    def __init__(self):
        file_name = Path(__file__).parent / "RHEA_MR_SBW_CPACS3.xml"
        super(CPACS_MR_SBW, self).__init__(file_name)

    def set_planform_design_variables(self):

        self.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                        "wing[@uID='MainWing']/sections/"
                                        "section[@uID='MainWing_sec1']/"
                                        "elements/element/transformation/"
                                        "scaling/x",
                                        name="root_chord",
                                        is_design_variable=True,
                                        lower_bound=None,
                                        upper_bound=None)

        self.select_twist_angle("twist_kink", "MainWing", "MainWing_sec2",
                                is_design_variable=True, 
                                lower_bound=None,
                                upper_bound=None)

        self.select_twist_angle("twist_tip", "MainWing", "MainWing_sec3",
                                is_design_variable=True,
                                lower_bound=None,
                                upper_bound=None)

        self.select_sweep_angle_variable("sweep_section_1",
                                         "MainWing",
                                         "MainWing_sec2_pos",
                                         is_design_variable=True,
                                         lower_bound=None,
                                         upper_bound=None)

    def set_planform_processed_variables(self):
        """
        Set variables that must be processed.

        Processed variables depends on external design variables.
        Processed variables must be processed when values are set into
        CPACS data.
        """

        # Tip chord
        tip_chord_proc = lambda root_chord, taper_ratio: taper_ratio*root_chord
        self.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                        "wing[@uID='MainWing']/sections/"
                                        "section[@uID='MainWing_sec3']/"
                                        "elements/element/transformation/"
                                        "scaling/x",
                                        "tip_chord",
                                        process=tip_chord_proc)

        # length section 1 and 2
        kink_location = 0.5
        section_length_1_proc = lambda span: kink_location*span
        self.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                        "wing[@uID='MainWing']/positionings/"
                                        "positioning[@uID='MainWing_sec2_pos']/length",
                                        "section_length_1",
                                        process=section_length_1_proc)

        section_length_2_proc = lambda span: span*(1-kink_location)
        self.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                        "wing[@uID='MainWing']/positionings/"
                                        "positioning[@uID='MainWing_sec3_pos']/length",
                                        "section_length_2",
                                        process=section_length_2_proc)

        # sweep section 2
        sweep_section_2_process = lambda sweep_section_1: sweep_section_1
        self.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                        "wing[@uID='MainWing']/positionings/"
                                        "positioning[@uID='MainWing_sec3_pos']/"
                                        "sweepAngle",
                                        "sweep_section_2",
                                        process=sweep_section_2_process)