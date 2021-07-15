

from typing import Mapping
import numpy as np
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).absolute().parents[2]))
from gemseo_rhea.cpacs import MDODisciplineCPACSBased
from gemseo_rhea.cpacs import CPACSStructureData

from RazorModel import RazorModel


class RazorDiscipline(MDODisciplineCPACSBased):
    """
    Definition of a RAZOR discipline for GEMSEO.

    The class inherits from MDODisciplineCPACSBased in order to be consitant
    with variables in CPACS file.
    In such way, all the inputs and the outputs that are used in the RAZOR
    discipline must exist in the CPACS file.
    """

    def __init__(self,
                 input_mapping,  # type: Mapping[str, xml.etree.ElementTree]
                 output_mapping,  # type: Mapping[str, xml.etree.ElementTree]
                 input_razor_template  # type: str
                 ):  # type: (...) -> None
        """
        Constructor.

        Args:
            input_mapping: The mapping to CPACS data of input variables.
            output_mapping: The mapping to CPACS data of output variables.
            input_razor_template: The name of the RAZOR input file template.
        """
        super(RazorDiscipline, self).__init__(input_cpacs_mapping=input_mapping,
                                              output_cpacs_mapping=output_mapping)

        self.__model = RazorModel(input_razor_template)

    def _run(self):  # type: (...) -> None
        """
        Run the model.

        This function overload the _run function of the GEMSEO MDODiscipline.
        It calls the model and store the results into local_data.
        """

        input_data = self.get_input_data()

        # since input_data associate variables to ndarray
        # we must convert to double (take the first value)
        input_data = {name: value[0] for name, value in input_data.items()}

        output_data = self.__model.execute(input_data)

        # again, output must be stored as ndarray
        for name, value in output_data.items():
            self.local_data[name] = np.array([value])


if __name__ == "__main__":

    # init the CPACS structure data
    cpacs_data = CPACSStructureData("razor_wing.xml")

    # select input
    cpacs_data.select_variable_from_xpath("./vehicles/aircraft/model/wings/wing/"
                                          "positionings/positioning[@uID='pos2']/"
                                          "sweepAngle",
                                          name="TIME")

    # select outputs
    cpacs_data.select_variable_from_xpath("./toolspecific/tool/Outputs/maxXMomentum",
                                          name="max_X_momentum")

    cpacs_data.select_variable_from_xpath("./toolspecific/tool/Outputs/"
                                          "meanYMomentum",
                                          name="mean_Y_momentum")

    # get the mappings that are used in the RAZOR discipline
    input_mapping = cpacs_data.get_sub_mapping(["TIME"])
    output_mapping = cpacs_data.get_sub_mapping(["max_X_momentum", "mean_Y_momentum"])

    input_file = str(Path(__file__).parent / "tutorial_naca0012_vol.cfg")
    
    # build the RAZOR discipline
    razor_discipline = RazorDiscipline(input_mapping,
                                       output_mapping,
                                       input_file)

    # GEMSEO uses ndarray even if the variable has only one dimension
    output = razor_discipline.execute(input_data={"TIME": np.array([0.01])})
    print(output)



