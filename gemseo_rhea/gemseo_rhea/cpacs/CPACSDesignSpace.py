

from gemseo.algos.design_space import DesignSpace


class CPACSDesignSpace(DesignSpace):
    """
    Build a design space based on CPACS data
    """

    def __init__(self, cpacs_mapping):
        super(CPACSDesignSpace, self).__init__()
        self._cpacs_mapping = cpacs_mapping

    def set_variable_data(self, name, l_b=None, u_b=None, value=None):
        """
        Add a new variable from the CPACS .xml file

        Args:
            name (str): variable name
            l_b (float or np.array): lower bound(s) (default is None)
            u_b (float or np.array): upper bound(s) (default is None)
            value (float or np.array): value(s) of variable (default is None,
                                       which means that the values are initialised
                                       from values read in .xml file)
        """
        try:
            size = self._cpacs_mapping.get_variable_size(name)
            if value is None:
                value = self._cpacs_mapping[name]
        except KeyError:
            raise KeyError("Variable {} does not exist in CPACS mapping".format(name))

        self.add_variable(name=name, size=size, l_b=l_b, u_b=u_b, value=value, var_type=DesignSpace.FLOAT)

    def set_cpacs_values(self):
        """
        Set the current design space values into cpacs struture.
        """
        for name_var in self.variables_names:
            self._cpacs_mapping[name_var] = self.get_current_x([name_var])
