
import numpy as np
import os
from typing import Mapping, Optional
import logging

from gemseo.utils.py23_compat import Path
from gemseo.api import configure_logger

import sharpy_model as sharpy

from gemseo_rhea.cpacs import MDODisciplineCPACSBased, CPACSStructureData

LOGGER = logging.getLogger()
configure_logger()


class SharpyModel:
    """Example of SHARPy model.

    This class is an example of the SHARPy model definition.
    """

    def __init__(self,
                 cpacs_file_name=None  # type: Optional[str]
                 ):  # type (...) -> None
        """Initialize the model.

        This function enables to initialize the model.
        For that purpose, an input CPACS file can be read
        in order to pick any variables up.

        Args:
            cpacs_file_name: the CPACS file
        """

        if cpacs_file_name is not None:
            # read cpacs file and store constant parameters values
            # ...
            pass
        else:
            LOGGER.info("SHARPy model: none CPACS file.")

    def execute(self,
                input_values  # type: Mapping[str, np.ndarray]
                ):  # type (...) -> Mapping[str, np.ndarray]
        """Execute the model.

        Args:
            input_values: input values

        Returns:
            A dictionary associating output names to output values
        """

        # get sweep value
        sweep = input_values["sweep"] * np.pi / 180

        # build the wing model based on sweep
        wing = sharpy.build_wing_data(sweep)

        # set sharpy solver and build input files
        sharpy.build_solvers_settings(wing)

        # get path and name of file
        file_path = Path(wing.case_route)
        file_name = Path(wing.case_name+".sharpy")

        # run sharpy
        abs_path_file = file_path / file_name
        os.system("sharpy {}".format(abs_path_file))

        # parse the output
        response = self.parse_results()

        return response

    def parse_results(self):  # type (...) -> Mapping[str, np.ndarray]
        """Parse any output file.

        Returns:
            A dictionary of responses
        """
        return {"mach_divergence": np.array([0.8])}


class SharpyDiscipline(MDODisciplineCPACSBased):
    """SHARPy discipline for GEMSEO.

    This class enables to build the discipline
    that is used inside MDO process.
    It uses a SharpyModel in order to init and
    run the model.
    """

    def __init__(self,
                 input_cpacs_mapping,  # type: CPACSMapping
                 output_cpacs_mapping,  # type: CPACSMapping
                 name_cpacs_file=None,  # type: Optional[str]
                 ):  # type (...) -> None
        """Initialization of the discipline.

        Args:
            input_cpacs_mapping: the input cpacs mapping
            output_cpacs_mapping: the output cpacs mapping
            name_cpacs_file: the name of cpacs file
        """

        # call the mother class constructor
        super(SharpyDiscipline, self).__init__(input_cpacs_mapping=input_cpacs_mapping,
                                               output_cpacs_mapping=output_cpacs_mapping,
                                               name="SHARPy_discipline")

        # build the shapry model
        self._sharpy_model = SharpyModel(name_cpacs_file)

    def _run(self):  # type (...) -> None
        """Run the discipline.

        This is called when executing discipline.
        It gets input values, runs the model,
        and stores the outputs to local_data attribute.
        """

        # get a dict of current values
        input_values = self.get_input_data()

        # compute the model
        output_values = self._sharpy_model.execute(input_values)

        # store outputs inside the discipline
        for name, val in output_values.items():
            self.local_data[name] = val

        return


if __name__ == "__main__":

    # select variables from CPACS
    cpacs_data = CPACSStructureData("sharpy_wing.xml")
    cpacs_data.select_variable_from_xpath(xpath="./vehicles/aircraft/model/wings"
                                             "/wing/positionings/"
                                             "positioning[@uID='pos3']/sweepAngle",
                                          name="sweep")

    cpacs_data.select_variable_from_xpath(xpath="./vehicles/aircraft/model/"
                                              "analyses/aeroelastics/"
                                              "divergence/cases/case/ma",
                                          name="mach_divergence")

    # build the discipline
    input_mapping = cpacs_data.get_sub_mapping(["sweep"])
    output_mapping = cpacs_data.get_sub_mapping(["mach_divergence"])
    sharpy_disc = SharpyDiscipline(input_mapping, output_mapping)

    # execute the discipline
    sharpy_disc.execute({"sweep": np.array([35])})

    # set current disciplines values to cpacs file and write a new file
    sharpy_disc.set_cpacs_values()
    cpacs_data.write_xml("out_sharpy_wing.xml")

