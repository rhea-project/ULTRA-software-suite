
import pytest
import numpy as np

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSStructureData, CPACSDesignSpace

dir_name = Path(__file__).parent.absolute()


def build_cpacs_struct():
    return CPACSStructureData(dir_name / "test_input_wing.xml")


def test_set_non_existing_design_space_var():
    """Test the addition of a non-existing var to design space"""
    cpacs_struct = build_cpacs_struct()
    cpacs_space = CPACSDesignSpace(cpacs_struct.input_mapping)
    with pytest.raises(KeyError, match="Variable my_var does not exist in CPACS mapping"):
        cpacs_space.set_variable_data("my_var", l_b=-10, u_b=10)


def test_add_one_cpacs_variable_size_1():
    """Test the addition of one variable from xpath."""
    cpacs_struct = build_cpacs_struct()
    name_var = "upper_naca"
    cpacs_struct.select_input_from_xpath("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB",
                                         name_var)

    cpacs_space = CPACSDesignSpace(cpacs_struct.input_mapping)
    cpacs_space.set_variable_data(name_var, l_b=-10, u_b=10)

    current_x = cpacs_space.get_current_x()
    variable_names = cpacs_space.get_indexed_variables_names()

    assert current_x == pytest.approx(0.2)
    assert name_var in variable_names


def test_add_two_cpacs_variable_size_1_3():
    """Test the addition of two variables from xpath with different sizes."""
    cpacs_struct = build_cpacs_struct()
    name_var_u = "upper_naca"
    cpacs_struct.select_input_from_xpath("./vehicles/profiles/"
                                         + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB", name_var_u)
    name_var_l = "lower_naca"
    cpacs_struct.select_input_from_xpath("./vehicles/profiles/"
                                          + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB", name_var_l)

    cpacs_space = CPACSDesignSpace(cpacs_struct.input_mapping)
    cpacs_space.set_variable_data(name_var_u)
    cpacs_space.set_variable_data(name_var_l)

    current_x = cpacs_space.get_current_x()
    variable_names = cpacs_space.get_indexed_variables_names()

    assert (current_x == pytest.approx(np.array([0.2, 0.15, 0.1, 0.001])))
    assert name_var_u in variable_names
    assert name_var_l + "!0" in variable_names
    assert name_var_l + "!1" in variable_names
    assert name_var_l + "!2" in variable_names
    assert len(variable_names) == 4


def test_set_cpacs_variable_size_1():
    """Test the addition of one variable from xpath with size 1."""
    cpacs_struct = build_cpacs_struct()
    name_var_u = "upper_naca"
    cpacs_struct.select_input_from_xpath("./vehicles/profiles/"
                                         + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB", name_var_u)

    cpacs_space = CPACSDesignSpace(cpacs_struct.input_mapping)
    cpacs_space.set_variable_data(name_var_u, value=999.)
    cpacs_space.set_cpacs_values()

    assert cpacs_space._cpacs_mapping._dict_elements[name_var_u].text == "999.0"


def test_set_cpacs_variable_size_3():
    """Test the addition of one variable from xpath with size 3."""
    cpacs_struct = build_cpacs_struct()
    name_var_l = "lower_naca"
    cpacs_struct.select_input_from_xpath("./vehicles/profiles/"
                                         + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                         name_var_l)

    cpacs_space = CPACSDesignSpace(cpacs_struct.input_mapping)
    cpacs_space.set_variable_data(name_var_l, value=np.array([11., 22., 33.]))

    cpacs_space.set_cpacs_values()

    assert cpacs_space._cpacs_mapping._dict_elements[name_var_l].text == "11.0;22.0;33.0"


