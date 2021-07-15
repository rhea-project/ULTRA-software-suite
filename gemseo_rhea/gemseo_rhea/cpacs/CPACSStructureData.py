

import inspect
import xml.etree.ElementTree as ET
from typing import Iterable

from .CPACSMapping import CPACSMapping
from .CPACSVariable import CPACSVariable


class CPACSStructureData:
    """

    """

    _PROFILE_CST_LOWER = CPACSVariable(
        "./vehicles/profiles/wingAirfoils/"
        "wingAirfoil[@uID='{airfoil_ID}']/cst2D/lowerB",
        ["airfoil_ID"])

    _PROFILE_CST_UPPER = CPACSVariable(
        "./vehicles/profiles/wingAirfoils/"
        "wingAirfoil[@uID='{airfoil_ID}']/cst2D/upperB",
        ["airfoil_ID"])

    _SWEEP_ANGLE = CPACSVariable(
        "./vehicles/aircraft/model/wings/wing[@uID='{wing_ID}']/"
        "positionings/positioning[@uID='{position_ID}']/sweepAngle",
        ["wing_ID", "position_ID"])

    _TWIST_ANGLE = CPACSVariable(
       "./vehicles/aircraft/model/wings/wing[@uID='{wing_ID}']/"
       "sections/section[@uID='{section_ID}']/transformation/rotation/y",
        ["wing_ID", "section_ID"])

    def __init__(self, file_name):
        self.__cpacs_mapping = CPACSMapping()
        self.__design_variables = set()
        self.__bounds = dict()
        self.__processed_variables = dict()
        self.__tree = ET.parse(file_name)

    @property
    def design_variables(self):  # type: (...) -> Iterable
        """Return the design variables."""
        return self.__design_variables

    @property
    def processed_variables(self):  # type: (...) -> Iterable
        """Return the processed variables."""
        return self.__processed_variables.keys()

    def get_process_arguments(self,
                              name  # type: str
                              ):  # type: (...) -> Iterable
        """
        Get the arguments of a process function.

        Args:
            name: The name of the processed variable.

        Returns:
            The arguments of the process function.
        """
        return inspect.signature(self.__processed_variables[name]).parameters.keys()

    def select_variable_from_xpath(self,
                                   xpath,  # type: str
                                   name,  # type: str
                                   is_design_variable=False,  # type:bool
                                   lower_bound=None,  # type: Optional[float]
                                   upper_bound=None,  # type: Optional[float]
                                   process=None  # type: Optional[callable]
                                   ):  # type (...) -> None
        """
        Select an input variable in xml tree from XPath

        Args:
            xpath: The xml Xpath of the selected variable.
            name: The name of the variable.
        """
        element = self.__get_xml_element(xpath)
        self.__cpacs_mapping.add_xml_element(name, element)

        self.__set_variable_properties(name, is_design_variable,
                                       lower_bound, upper_bound,
                                       process)

    def __set_variable_properties(self, name,
                                  is_design_variable=False,
                                  lower_bound=None,
                                  upper_bound=None,
                                  process=None  # type: Optional[callable]
                                  ):  # type: (...) -> None
        """
        Set properties of design variables.

        Args:
            process:
            name: The name of the variable.
            is_design_variable: True is the variable is a design variable.
            lower_bound: The lower bound value (only for design variables).
            upper_bound: The upper bound value (only for design variables).
            process: The process executed when the variable is set with new values.
        """
        if is_design_variable:
            self.__design_variables.add(name)
            self.__bounds.update({name: (lower_bound, upper_bound)})

        if process:
            if is_design_variable:
                raise ValueError("CPACS variable cannot be a design variable"
                                 "and a processed variable in the same time.")

            self.__processed_variables.update({name: process})

    def __get_xml_element(self, xpath):
        """
        Return the xml element corresponding to prescribed XPath

        Args:
            xpath:

        Returns:

        """
        elements = self.__xml_find_xpath(xpath)
        # check that xpath contains exactly 1 element
        if len(elements) == 0:
            raise ValueError("None element found corresponding "
                             + "to the current XPath.")
        elif len(elements) > 1:
            raise ValueError("Current XPath has more than "
                             + "one element.")
        else:
            element = elements[0]

        return element

    def __xml_find_xpath(self, xpath):
        """
        Find XPath in xml tree

        Args:
            xpath (str): the required XPath

        Returns:
            A list of xml elements
        """
        return self.__tree.findall(xpath)

    def write_xml(self, name_file):
        """
        Write .xml file from data

        Args:
            name_file (str): the file name
        """
        self.__tree.write(name_file)

    def __str__(self):
        """
        Convert object to string

        Returns:

        """
        return "CPACS VARIABLES =\n{}\n".format(str(self.__cpacs_mapping))

    def get_value(self, name_var):
        """
        Return the values corresponding to the required input variable

        Args:
            name_var: The name of variable.

        Returns:

        """
        return self.__cpacs_mapping[name_var]

    def set_value(self,
                  name_var,  # type: str
                  **kwargs  # type: float
                  ):  # type: (...) -> None
        """
        Set new values to the variable.

        Args:
            name_var: The name of the variable.
            kwargs: The values of arguments.
                If the variable is a processed variable, kwargs must include
                all keyword arguments needed in the processed function.
                Otherwise, only keyword argument ``value`` should be provided.
        """
        process = self.__processed_variables.get(name_var)
        if process:
            # here we build a sub_data which include only the values
            # needed in the process function, because kwargs could
            # include more useless data.
            # process_args = inspect.signature(process).parameters.keys()
            # sub_data = {var: kwargs[var] for var in process_args}
            self.__cpacs_mapping[name_var] = process(**kwargs)
        else:
            self.__cpacs_mapping[name_var] = kwargs["value"]

    def get_xml_element(self,
                        name  # type: str
                        ):  # type (...) -> xml.etree.ElementTree
        """Get the xml element of the variable.

        Args:
            name: The name of the variable.

        Returns:
            The xml element tree.
        """
        return self.__cpacs_mapping.get_xml_element(name)

    def get_variable_size(self,
                          name  # type: str
                          ):  # type: (...) -> int
        """Get the size of the variable.

        Args:
            name: The nme of the variable

        Returns:
            The size of the variable
        """
        return self.__cpacs_mapping.get_variable_size(name)

    def get_bounds(self,
                   name  # type: str
                   ): # type: (...) -> Tuple[float, float]
        """Get the bounds of the design variable.

        Selected variables are design variables only if
        is_design_variable=True.

        Args:
            name: The name of the design variable.

        Returns:
            The values of bounds.
        """
        return self.__bounds[name]

    def get_sub_mapping(self,
                       variables  # type: Iterable[str]
                       ):  # type (...) -> Mapping[str, ET]
        """
        Get a sub-mapping.

        Args:
            variables:

        Returns:

        """
        sub_mapping = CPACSMapping()

        for variable in variables:
            sub_mapping.add_xml_element(variable,
                                        self.__cpacs_mapping.get_xml_element(variable))

        return sub_mapping

    def __iter__(self):
        """Iterate over the selected design variables with its value.

        This function enables to iterate over the set
        of design variables (only, not all selected variables),
        returning its name and value.
        Design variables are variables selected with is_design_variable=True.

        Examples:
            >>> cpacs = CPACSStructureData(file_name)
            >>> # select any variables...
            >>> for var, val in cpacs:
            >>>    print(var, val)
        """
        for var in self.__design_variables:
            yield var, self.get_value(var)

    def select_twist_angle(self,
                           name,
                           wing_ID,
                           section_ID,
                           is_design_variable=False,
                           lower_bound=None,
                           upper_bound=None):
        """Select the twist angle.

        Args:
            name: The name of the variable.
            wing_ID: The wing ID.
            section_ID: The section ID.
        """
        xpath = self._TWIST_ANGLE.get_xpath([wing_ID, section_ID])
        self.select_variable_from_xpath(xpath, name, is_design_variable,
                                        lower_bound, upper_bound)

    def select_lower_profile_variable(self,
                                      name,
                                      airfoil_ID,
                                      is_design_variable=False,  # type:bool
                                      lower_bound=None,  # type: float
                                      upper_bound=None,  # type: float
                                      ):
        """
        Add variable that corresponds to lower profile airfoil.

        Args:
            name (str): variable name
            airfoil_ID (str): the name of airfoil ID (@uID xml attribute)
        """
        xpath = self._PROFILE_CST_LOWER.get_xpath([airfoil_ID])
        self.select_variable_from_xpath(xpath, name,
                                        is_design_variable, lower_bound, upper_bound)

    def select_upper_profile_variable(self,
                                      name,
                                      airfoil_ID,
                                      is_design_variable=False,  # type:bool
                                      lower_bound=None,  # type: float
                                      upper_bound=None,  # type: float
                                      ):
        """
        Add variable that corresponds to upper profile airfoil

        Args:
            name (str): variable name
            airfoil_ID (str): the name of airfoil ID (@uID xml attribute)
        """
        xpath = self._PROFILE_CST_UPPER.get_xpath([airfoil_ID])
        self.select_variable_from_xpath(xpath, name,
                                        is_design_variable, lower_bound, upper_bound)

    def select_sweep_angle_variable(self,
                                    name,
                                    wing_ID,
                                    position_ID,
                                    is_design_variable=False,  # type:bool
                                    lower_bound=None,  # type: float
                                    upper_bound=None,  # type: float
                                    ):
        """
        Add variable that corresponds to wing section sweep angle

        Args:
            name (str): variable name
            wing_ID (str): the name of a wing ID (@uID xml attribute)
            position_ID (str): the name of a position ID (@uID xml attribute)
        """
        xpath = self._SWEEP_ANGLE.get_xpath([wing_ID, position_ID])
        self.select_variable_from_xpath(xpath, name, is_design_variable,
                                        lower_bound, upper_bound)




