

from copy import deepcopy, copy
from typing import Optional

from gemseo.algos.design_space import DesignSpace
from gemseo_rhea.cpacs.CPACSStructureData import CPACSStructureData
from gemseo_rhea.cpacs.CPACSMapping import CPACSMapping


class CPACSDesignSpace(DesignSpace):
    """
    Build a design space based on CPACS data.
    """

    def __init__(self,
                 cpacs_structure_data=None  # type: Optional[CPACSStructureData]
                 ):
        """

        Args:
            cpacs_structure_data: The CPACS structure data.
                If ``None``, the design space will be empty (no
                design variable).
        """
        super(CPACSDesignSpace, self).__init__()

        self._cpacs_structure = cpacs_structure_data

        self.__set_design_variable()

    def __set_design_variable(self):  # type: (...) -> None
        """
        Set variable to the deisgn space
        """
        for variable, value in self._cpacs_structure:
            size = self._cpacs_structure.get_variable_size(variable)
            bounds = self._cpacs_structure.get_bounds(variable)

            self.add_variable(variable, size=size,
                              value=value, l_b=bounds[0], u_b=bounds[1])

    def __deepcopy__(self, memo={}):
        """This function is called when doing deepcopy.

        This function is overridden because
        we don't want __cpacs_mapping to be deepcopied
        when doing deepcopy of the whole object.
        Indeed, __cpacs_mapping is a reference toward the cpacs structure
        so we want to keep that reference when executing :func:`set_cpacs_value`.
        This means that if a design space is deepcopied,
        executing :func:`set_cpacs_value` on one after the other will override
        the data.
        Design space deepcopy should be used in order to filter design variables
        such as two spaces should involve different variables.

        Returns:
            A new object instance.
        """
        new_obj = type(self)(self._cpacs_structure)

        for attr_name, attr_value in self.__dict__.items():
            if attr_name is not "_cpacs_structure":
                new_obj.__dict__[attr_name] = deepcopy(attr_value)

        return new_obj

    def set_cpacs_values(self):
        """
        Set the current design space values into cpacs structure.
        """
        for var in self.variables_names:
            if var in self._cpacs_structure.design_variables:
                self._cpacs_structure.set_value(var, value=self.get_current_x([var]))

        data = self.get_current_x_dict()
        for var in self._cpacs_structure.processed_variables:
            args = self._cpacs_structure.get_process_arguments(var)

            # we make sure that all args are included in the current design space,
            # otherwise process arguments have been splitted.
            sub_data = {arg: data[arg] for arg in args if arg in data.keys()}

            if len(sub_data) != len(args):
                raise ValueError("Some arguments of the function related to the "
                                 "processed variable {} "
                                 "are missing in the current design space. "
                                 "Needed arguments are {}.".format(var, list(args)))

            self._cpacs_structure.set_value(var, **sub_data)


        # for var in self.__cpacs_structure.design_variables:
        #     self.__cpacs_structure.set_value(var, value=self.get_current_x([var]))
        # for var in self.__cpacs_structure.processed_variables:
        #     self.__cpacs_structure.set_value(var, self.get_current_x_dict())
