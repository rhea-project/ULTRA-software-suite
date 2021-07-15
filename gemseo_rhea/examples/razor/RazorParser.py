
from enum import Enum
import numpy as np
import h5py
from pathlib import Path

class RazorParser:
    """
    Parser for the RAZOR input file.
    """

    class TargetFormat(Enum):
        """
        Format of target output file.
        """
        SU2 = "SU2_CSV"
        PARAVIEW = "PARAVIEW"
        TECLPLOT = "TECPLOT"

    def __init__(self,
                 input_file  # type: str
                 ):  # type: (...) -> None
        """
        Constructor.

        Args:
            input_file: The name of the input file.
        """
        self.__input_file = input_file

        # name of database
        self.__database = ""
        # name of parameters in target list
        self.__variables = []
        # name of fields of interest for output
        self.__fields = []
        # format of output results
        self.__target_format = RazorParser.TargetFormat.SU2

        self.__parse_input()

        self.__check()

    def __check(self):  # type: (...) -> None
        """
        Check that all attributes are correctly set.
        """
        if self.__database == "":
            raise ValueError("Database has bot been found.")
        if len(self.__variables) == 0:
            raise ValueError("None variable have been found.")
        if len(self.__fields) == 0:
            raise ValueError("None fields have been found.")

    @property
    def input_file(self):
        return self.__input_file
        
    @property
    def variables(self):  # type: (...) -> List[str]
        """Return variables."""
        return self.__variables

    @property
    def fields(self):  # type: (...) -> List[str]
        """Return fields."""
        return self.__fields

    def is_training_in_database(self):  # type: (...) -> bool
        """
        Check if training data are in the database.

        Returns:
            True if training data exists, False otherwise.
        """
        if not Path(self.__database).exists():
            return False
        
        f = h5py.File(self.__database, "r")
        train = f["TrainingSet"]
        if len(train) > 0:
            f.close()
            return True
        else:
            f.close()
            return False

    def is_model_in_database(self):  # type: (...) -> None
        """
        Check if the model is in database.

        Returns:
            True if the model exist, False otherwise.
        """
        if not Path(self.__database).exists():
            return False

        f = h5py.File(self.__database, "r")
        model = f["LDmodels"]
        if len(model) > 0:
            f.close()
            return True
        else:
            f.close()
            return False

    def __parse_input(self):  # type: (...) -> None
        """
        Parse the input file.

        This function fills attributes.
        """
        f = open(self.__input_file, "r")

        for line in f:
            if self.__get_tag(line) == "NAME_PARM":
                self.__variables = self.__get_values(line)
            if self.__get_tag(line) == "FIELDS_NAME":
                self.__fields = self.__get_values(line)
            if self.__get_tag(line) == "DATABASE_NAME":
                self.__database = self.__get_values(line)[0]

        f.close()

    def __get_tag(self,
                  text  # type: str
                  ):  # type: (...) -> str
        """
        Get the tag of a line before ``=`` sign.

        Args:
            text: The line of the file.

        Returns:
            The tag.
        """
        return text.split("=")[0].strip()

    def __get_values(self,
                     text  # type: str
                     ):  # type: (...) -> List[str]
        """
        Get the values found in a line.

        Args:
            text: The line of the file.

        Returns:
            All values found.
        """
        values = text.split("=")[1].split(",")
        return [value.strip() for value in values]

    def __build_target_list(self,
                            values,  # type: Mapping[str, double]
                            file,  # type: str
                            ):  # type: (...) -> str
        """
        Build the target list line.

        Args:
            values: The values of parameters of the target list.
            file: The file where predicted data are exported.

        Returns:
            The text line of the target list.
        """
        target_txt = "TARGET_LIST="
        try:
            values_txt = ",".join([str(values[var]) for var in self.variables])
        except KeyError:
            raise NameError("Some expected variables are not found in provided "
                            "variables.\n"
                            "\t- Expected variables are {}.\n"
                            "\t- Provides variables are {}.\n".format(self.variables,
                                                                    list(values.keys())))
        return target_txt + values_txt + "," + file+"\n"

    def generate_target_file(self,
                             target_values,  # type: Mapping[str, double]
                             output_target_file,  # type: str
                             generated_file,  # type: str
                             target_format=TargetFormat.SU2  # type: RazorParser.TargetFormat
                             ):  # type: (...) -> None
        """
        Generate a new input file with specified target list values.

        Args:
            target_values: The target list of values.
            output_target_file: The file were predicted data are exported.
            generated_file: The name of the new generated input file.
            target_format: The format of the file for the exported data.
        """

        self.__target_format = target_format

        new_f = open(generated_file, "w")
        f = open(self.__input_file, "r")

        for line in f:
            if "TARGET_LIST=" in line:
                new_f.write(self.__build_target_list(target_values, output_target_file))
            elif "TARGET_FMT=" in line:
                new_f.write("TARGET_FMT={}\n".format(target_format.value))
            else:
                new_f.write(line)

        f.close()
        new_f.close()

    def parse_output(self,
                     name_file  # type: str
                     ):  # type: (...) -> Mapping[str, ndarray]
        """
        Parse the output file containing targeted results.

        Args:
            name_file: The name of the output file.

        Returns:
            The values of the output file.
        """

        if self.__target_format is not RazorParser.TargetFormat.SU2:
            raise ValueError("Output format {} is not implemented "
                             "for parsing.".format(RazorParser.TargetFormat.SU2.value))

        nb_col = len(self.fields)

        f = open(name_file, "r")
        header = f.readline().split(",")[0:nb_col]
        data = np.loadtxt(f, delimiter=",", usecols=range(nb_col))
        f.close()

        header = [self.__clean_header_name(name) for name in header]

        # check that header correspond to required fields
        for i, field in enumerate(self.fields):
            if field != header[i]:
                raise ValueError("Field {} is not found in output file.".format(field))

        return {name: data[:, i] for i, name in enumerate(header)}

    def __clean_header_name(self,
                            name  # type: str
                            ):  # type: (...) -> str
        """
        Clean the name of output fields in header.

        Args:
            name: The name of the field.

        Returns:
            The new cleaned name.

        """
        name = name.replace('"', '')
        name = name.replace("'", "")
        return name

        
