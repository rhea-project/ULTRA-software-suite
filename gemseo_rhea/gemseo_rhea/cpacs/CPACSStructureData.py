
import xml.etree.ElementTree as ET

from .CPACSMapping import CPACSMapping
from .CPACSVariable import CPACSVariable


class CPACSStructureData:
    """

    """

    PROFILE_CST_LOWER = CPACSVariable(
        "./vehicles/profiles/wingAirfoils/"
        +"wingAirfoil[@uID='{airfoil_ID}']/cst2D/lowerB",
        ["airfoil_ID"])

    PROFILE_CST_UPPER = CPACSVariable(
        "./vehicles/profiles/wingAirfoils/"
        +"wingAirfoil[@uID='{airfoil_ID}']/cst2D/upperB",
        ["airfoil_ID"])

    SWEEP_ANGLE = CPACSVariable(
        "./vehicles/aircraft/model/wings/wing[@uID='{wing_ID}']/"
        +"positionings/positioning[@uID='{position_ID}']/sweepAngle",
        ["wing_ID", "position_ID"])

    def __init__(self, file_name):
        self._input_mapping = CPACSMapping()
        self._output_mapping = CPACSMapping()
        self._tree = ET.parse(file_name)

    def select_input_from_xpath(self, xpath, name):
        """
        Select an input variable in xml tree from XPath

        Args:
            xpath:
            name:

        Returns:

        """
        element = self._get_element(xpath)
        self._input_mapping.add_xml_element(name, element)

    def select_output_from_xpath(self, xpath, name):
        """
        Select an output response in xml tree from XPath

        Args:
            xpath:
            name:

        Returns:

        """
        element = self._get_element(xpath)
        self._output_mapping.add_xml_element(name, element)

    def _get_element(self, xpath):
        """
        Return the xml element corresponding to prescribed XPath

        Args:
            xpath:

        Returns:

        """
        elements = self._xml_find_xpath(xpath)
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

    def _xml_find_xpath(self, xpath):
        """
        Find XPath in xml tree

        Args:
            xpath (str): the required XPath

        Returns:
            A list of xml elements
        """
        return self._tree.findall(xpath)

    def write_xml(self, name_file):
        """
        Write .xml file from data

        Args:
            name_file (str): the file name
        """
        self._tree.write(name_file)

    def __str__(self):
        """
        Convert object to string

        Returns:

        """
        str = ""
        str += "* INPUTS =\n"+str(self._input_mapping)+"\n"
        str += "* OUTPUTS =\n"+str(self._output_mapping)
        return str

    def get_input_values(self, name_var):
        """
        Return the values corresponding to the required input variable

        Args:
            name_var:

        Returns:

        """
        return self._input_mapping[name_var]

    def set_input_values(self, name_var, values):
        """
        Set new values for a prescribed input variable

        Args:
            name_var:
            values:

        Returns:

        """
        self._input_mapping[name_var] = values

    def get_output_values(self, name_var):
        """
        Return the values of a required output variable

        Args:
            name_var:

        Returns:

        """
        return self._ouput_mapping[name_var]

    def set_output_values(self, name_var, values):
        """
        Set new values for the prescribed output variable

        Args:
            name_var:
            values:

        Returns:

        """
        self._output_mapping[name_var] = values

    @property
    def input_mapping(self):
        return self._input_mapping

    @property
    def output_mapping(self):
        return self._output_mapping

    def select_lower_profile_variable(self, name, airfoil_ID):
        """
        Add variable that corresponds to lower profile airfoil

        Args:
            name (str): variable name
            airfoil_ID (str): the name of airfoil ID (@uID xml attribute)
        """
        xpath = self.PROFILE_CST_LOWER.get_xpath([airfoil_ID])
        self.select_input_from_xpath(xpath, name)

    def select_upper_profile_variable(self, name, airfoil_ID):
        """
        Add variable that corresponds to upper profile airfoil

        Args:
            name (str): variable name
            airfoil_ID (str): the name of airfoil ID (@uID xml attribute)
        """
        xpath = self.PROFILE_CST_UPPER.get_xpath([airfoil_ID])
        self.select_input_from_xpath(xpath, name)

    def select_sweep_angle_variable(self, name, wing_ID, position_ID):
        """
        Add variable that corresponds to wing section sweep angle

        Args:
            name (str): variable name
            wing_ID (str): the name of a wing ID (@uID xml attribute)
            position_ID (str): the name of a position ID (@uID xml attribute)
        """
        xpath = self.SWEEP_ANGLE.get_xpath([wing_ID, position_ID])
        self.select_input_from_xpath(xpath, name)