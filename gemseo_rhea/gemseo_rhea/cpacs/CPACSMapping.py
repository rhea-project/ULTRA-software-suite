
import numpy as np
import xml.etree.ElementTree


class CPACSMapping:
    """
    CPACS Mapping

    """

    def __init__(self):
        # associate variables to xml tree elements (enables to set xml)
        self.__dict_elements = {}

    def _xml_get_text(self, name_var):
        """
        Get the text corresponding to a required already stored variable

        Args:
            name_var (str): the variable name

        Returns:
            The string element text
        """
        return self.__dict_elements[name_var].text

    def _xml_set_text(self, name_var, text_value):
        """
        Set the text corresponding to a required already stored variable

        Args:
            name_var (str): the variable name
            text_value (str): the text
        """
        self.__dict_elements[name_var].text = text_value

    def add_xml_element(self,
                        var_name,  # type: str
                        element  # type: xml.etree.ElementTree
                        ):  # type: (...) -> None
        """Add a xml element.

        Args:
            var_name: The name of the variable.
            element: The xml element.
        """
        self.__dict_elements.update({var_name : element})

    def get_xml_element(self,
                        var_name  # type: str
                        ):  # type: (...) -> xml.etree.ElementTree
        """Return a xml element.

        Args:
            var_name: The nmae of the variable

        Returns:
            An xml element.
        """
        return self.__dict_elements[var_name]

    def __getitem__(self,
                    name_var  # type: str
                    ):  # type (...) -> np.ndarray
        """
        Get the value corresponding to the variable.

        Args:
            name_var: The variable name.

        Returns:
            The value of the variable.
        """
        element_text = self._xml_get_text(name_var)
        values = np.array([float(x) for x in element_text.split(";")])
        return values

    def __setitem__(self,
                    name_var,  # type: str
                    values  # type: Tuple[float, np.ndarray]
                    ):  # type (...) -> None
        """
        Set new value to the variable.

        Args:
            name_var: The name of the variable.
            values: The values.
        """
        if type(values) is float or type(values) is int:
            values = np.array([values])

        if self.get_variable_size(name_var) != len(values):
            raise ValueError("Variable {}. Prescribed "
                             "values has the wrong size.".format(name_var))

        list_values = [str(x) for x in values]
        text_values = ";".join(list_values)
        self._xml_set_text(name_var, text_values)

    def get_variable_size(self, name_var):
        """
        Check the variable size

        Args:
            name_var (str): the variable name

        Returns:
            integer: the size of variable
        """
        return len(self.__getitem__(name_var))

    def __str__(self):
        """
        String representation

        Returns:
            string: the string representation
        """
        loc_str = ""
        for k, v in self.__dict_elements.items():
            loc_str += "- {}\n\t{}\n".format(k, v.text)

        return loc_str

    def get_dict_with_values(self):
        """"
        Return a dictionary of variables mapped with array of values
        """
        return {name: self[name] for name in self.__dict_elements.keys()}

    def keys(self):  # type: (...) -> List
        """Get mapping keys."""
        return list(self.__dict_elements.keys())