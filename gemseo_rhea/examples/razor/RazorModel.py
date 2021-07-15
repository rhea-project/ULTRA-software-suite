
import os
import numpy as np
from typing import Mapping
from enum import Enum
import sys

from gemseo.utils.py23_compat import Path

sys.path.append(str(Path(__file__).absolute().parent))
from RazorParser import RazorParser


class RazorModel:
    """
    This class define the RAZOR model.

    It is an example of how RAZOR could be interfaced in Python.

    Examples:
        >>> razor = RazorModel("input_template.cfg")
        >>> output = razor.execute({"parameter_1": 0.1, "parameter_2": 10.2})
    """

    # Environment variable for the RAZOR directory
    __RAZOR_ENV = "RAZOR_ENV"

    class RazorOptions(Enum):
        """Razor execution options."""
        PRE_PROC = "-p"
        BUILD = "-g"
        PREDICT = "-s"

    def __init__(self,
                 input_file,  # type: str
                 generated_target_file="generated_input_file.cfg",
                 output_file_prediction="output_results.csv"  # type: str
                 ):  # type: (...) -> None
        """
        Constructor.

        Args:
            input_file: The name of the input file.

        Raises:
            ValueError: If the environment variable of RAZOR directory is not set.
        """
        exec_path = os.getenv(self.__RAZOR_ENV)
        if exec_path is None:
            raise ValueError("Environment variable {} is not "
                             "set.".format(self.__RAZOR_ENV))

        self.__path_to_exec = Path(exec_path)
        self.__parser = RazorParser(input_file)
        self.__generated_target_file = generated_target_file
        self.__output_file_prediction = output_file_prediction

        # run preproc and model building only if they are not found
        # in database

        if not self.__parser.is_training_in_database():
            self.__run_razor(RazorModel.RazorOptions.PRE_PROC)
        if not self.__parser.is_model_in_database():
            self.__run_razor(RazorModel.RazorOptions.BUILD)

    def __post_treat_results(self,
                             values  # type: Mapping[str, ndarray]
                             ):  # type: (...) -> Mapping[str, double]
        """
        Post-treatment of outputs.

        This function enables to do any post-treatment on fields that
        are returned by RAZOR.

        Args:
            values: The values that must be post-treated.

        Returns:
            The new output values.
        """

        max_X = np.max(values["X-Momentum"])
        mean_Y = np.mean(values["Y-Momentum"])

        return {"max_X_momentum": max_X,
                "mean_Y_momentum": mean_Y}

    def __run_razor(self,
                    option  # type: RazorModel.RazorOptions
                    ):  # type: (...) -> None
        """
        Run the RAZOR model.

        Call the executable of RAZOR.

        Raises:
            ValueError: If the input file that must be run is not found.
        """
        executable = str(self.__path_to_exec / "razor")

        # the file that must executed depends on option...
        # for PREDICT case it must be the new generated file,
        # whereas for other cases it the basic inut_file
        if option == RazorModel.RazorOptions.PRE_PROC:
            input_file = self.__parser.input_file
        elif option == RazorModel.RazorOptions.BUILD:
            input_file = self.__parser.input_file
        elif option == RazorModel.RazorOptions.PREDICT:
            input_file = self.__generated_target_file
        else:
            raise ValueError("Option is unknown.")
        
        cmd = " ".join([executable, option.value, input_file])

        if not Path(input_file).exists():
            raise ValueError("The generated input file {} does "
                             "not exist.".format(input_file))

        os.system(cmd)

    def execute(self,
                values  # type: Mapping[str, double]
                ):  # type: (...) -> Mapping[str, ndarray]
        """
        Execute the offline model.

        This is the main call for the model.

        Args:
            values: The input values.

        Returns:
            The output values.
        """
        self.__parser.generate_target_file(values,
                                           self.__output_file_prediction,
                                           self.__generated_target_file)

        self.__run_razor(RazorModel.RazorOptions.PREDICT)

        output = self.__parser.parse_output(self.__output_file_prediction)

        post_treated_output = self.__post_treat_results(output)

        return post_treated_output

        

if __name__ == "__main__":

    input_file = "tutorial_naca0012_vol.cfg"

    razor_model = RazorModel(input_file)

    output = razor_model.execute({"TIME": 0.010})
    print(output)
