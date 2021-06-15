
import numpy as np


class CPACSMapping:
    """
    CPACS Mapping

    """

    def __init__(self):
        # associate variables to xml tree elements (enables to set xml)
        self._dict_elements = {}

    def _xml_get_text(self, name_var):
        """
        Get the text corresponding to a required already stored variable

        Args:
            name_var (str): the variable name

        Returns:
            The string element text
        """
        return self._dict_elements[name_var].text

    def _xml_set_text(self, name_var, text_value):
        """
        Set the text corresponding to a required already stored variable

        Args:
            name_var (str): the variable name
            text_value (str): the text
        """
        self._dict_elements[name_var].text = text_value

    def add_xml_element(self, var_name, element):
        self._dict_elements.update({var_name : element})

    def __getitem__(self, name_var):
        """
        Get values corresponding to the prescribed variable

        Args:
            name_var (str): the variable name

        Returns:
            numpy.array (float): required values
        """
        element_text = self._xml_get_text(name_var)
        values = np.array([float(x) for x in element_text.split(";")])
        return values

    def __setitem__(self, name_var, values):
        """
        Set new values to variable

        Args:
            name_var (str): the variable name
            values (numpy.array): array of values
        """
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
        for k, v in self._dict_elements.items():
            loc_str += "- "+ k + "\n" + "\t" + v.text + "\n"

        return loc_str

    def get_dict_with_values(self):
        """"
        Return a dictionnary of variables mapped with array of values
        """
        return {name: self[name] for name in self._dict_elements.keys()}

