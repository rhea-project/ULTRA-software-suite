
from gemseo_rhea.cpacs.CPACSDesignSpace import CPACSDesignSpace
from gemseo_rhea.elementary_bilevel.cpacs import CPACS_MR_SBW


class DesignSpace_MR_SBW(CPACSDesignSpace):

    def __init__(self, cpacs_struct):
        super(DesignSpace_MR_SBW, self).__init__(cpacs_struct)

    def add_planform_variables(self):

        initial_span = self._cpacs_structure.get_value("section_length_1") + \
                       self._cpacs_structure.get_value("section_length_2")

        self.add_variable("span", value=initial_span,
                          l_b=None, u_b=None)

        initial_taper_ratio = self._cpacs_structure.get_value("tip_chord") / \
                              self._cpacs_structure.get_value("root_chord")

        self.add_variable("taper_ratio", value=initial_taper_ratio,
                          l_b=None, u_b=None)

