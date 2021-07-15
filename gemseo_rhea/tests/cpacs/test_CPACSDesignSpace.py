
import pytest
import numpy as np
from copy import deepcopy

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSStructureData, CPACSDesignSpace

dir_name = Path(__file__).parent.absolute()


@pytest.fixture()
def cpacs_struct():
    return CPACSStructureData(dir_name / "test_input_wing.xml")


def test_add_one_cpacs_variable_size_1(cpacs_struct):
    """Test the addition of one variable from xpath."""
    name_var = "upper_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/"
                                            "wingAirfoils/wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB",
                                            name_var, is_design_variable=True,
                                            lower_bound=-10, upper_bound=10)

    cpacs_space = CPACSDesignSpace(cpacs_struct)

    current_x = cpacs_space.get_current_x()
    variable_names = cpacs_space.get_indexed_variables_names()
    lb = cpacs_space.get_lower_bound(name_var)
    ub = cpacs_space.get_upper_bound(name_var)

    assert current_x == pytest.approx(0.2)
    assert name_var in variable_names
    assert lb == pytest.approx(-10)
    assert ub == pytest.approx(10)


def test_setting_value_after_definition(cpacs_struct):
    """Test that changing design space value after definition is
    correctly taken into account into cpacs."""
    name_var = "upper_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/wingAirfoils/"
                                            "wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB",
                                            name_var, is_design_variable=True,
                                            lower_bound=-10,
                                            upper_bound=10)

    cpacs_space = CPACSDesignSpace(cpacs_struct)
    cpacs_space.set_current_variable(name_var, np.array([0.4]))

    cpacs_space.set_cpacs_values()

    assert cpacs_struct.get_value(name_var) == pytest.approx(0.4)


def test_add_two_cpacs_variable_size_1_3(cpacs_struct):
    """Test the addition of two variables from xpath with different sizes."""
    name_var_u = "upper_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/wingAirfoils/"
                                            "wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB", name_var_u,
                                            is_design_variable=True)
    name_var_l = "lower_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/wingAirfoils/"
                                            "wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/lowerB", name_var_l,
                                            is_design_variable=True)

    cpacs_space = CPACSDesignSpace(cpacs_struct)

    current_x_u = cpacs_space.get_current_x([name_var_u])
    current_x_l = cpacs_space.get_current_x([name_var_l])
    variable_names = cpacs_space.get_indexed_variables_names()

    assert current_x_l == pytest.approx(np.array([0.15, 0.1, 0.001]))
    assert current_x_u == pytest.approx(np.array([0.2]))
    assert name_var_u in variable_names
    assert name_var_l + "!0" in variable_names
    assert name_var_l + "!1" in variable_names
    assert name_var_l + "!2" in variable_names
    assert len(variable_names) == 4


def test_set_cpacs_variable_size_1(cpacs_struct):
    """Test the addition of one variable from xpath with size 1."""
    name_var_u = "upper_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/wingAirfoils/"
                                            "wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB", name_var_u,
                                            is_design_variable=True)

    cpacs_space = CPACSDesignSpace(cpacs_struct)
    cpacs_space.set_current_variable(name_var_u, np.array([999.]))

    cpacs_space.set_cpacs_values()

    assert cpacs_struct.get_value(name_var_u) == pytest.approx(999.)


def test_set_cpacs_variable_size_3(cpacs_struct):
    """Test the addition of one variable from xpath with size 3."""
    name_var_l = "lower_naca"
    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/"
                                            "wingAirfoils/wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/lowerB",
                                            name_var_l, is_design_variable=True)

    cpacs_space = CPACSDesignSpace(cpacs_struct)
    cpacs_space.set_current_variable(name_var_l, np.array([11., 22., 33.]))

    cpacs_space.set_cpacs_values()

    assert cpacs_struct.get_value(name_var_l) == pytest.approx(
        np.array([11.0, 22.0, 33.0]))


def test_design_space_deep_copy(cpacs_struct):
    """Test design space deepcopy.

    Make sure that when deepcopying, the internal cpacs mapping
    is still a reference.
    """
    name_var = "rib_rot"
    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/"
                                         "model/wings/wing/componentSegments/"
                                         "componentSegment/structure/"
                                         "ribsDefinitions/ribsDefinition/"
                                         "ribsPositioning/ribRotation/z",
                                            name_var, is_design_variable=True)
    init_cpacs_value = cpacs_struct.get_value(name_var)

    cpacs_main_space = CPACSDesignSpace(cpacs_struct)
    cpacs_main_space.set_current_variable(name_var, np.array([30.]))

    cpacs_copy_space = deepcopy(cpacs_main_space)
    cpacs_copy_space.set_current_variable(name_var, np.array([45.]))
    cpacs_copy_space.set_cpacs_values()

    assert init_cpacs_value == pytest.approx(90.)
    assert cpacs_main_space.get_current_x([name_var]) == pytest.approx(30.)
    assert cpacs_copy_space.get_current_x([name_var]) == pytest.approx(45.)
    assert cpacs_struct.get_value(name_var) == pytest.approx(45.)


def test_filter_on_design_space_deep_copy(cpacs_struct):
    """Test filtering on design space deepcopy.

    We make sure here that,
    if a deepcopied design space is filtered,
    changing values of variables, then setting cpacs values,
    is correctly handled.
    """
    name_var_1 = "x1"
    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/"
                                         "model/wings/wing/componentSegments/"
                                         "componentSegment/structure/"
                                         "ribsDefinitions/ribsDefinition/"
                                         "ribsPositioning/ribRotation/z",
                                            name_var_1, is_design_variable=True)

    name_var_2 = "x2"
    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                         "wing/positionings/positioning[3]/sweepAngle",
                                            name_var_2, is_design_variable=True)

    cpacs_main_space = CPACSDesignSpace(cpacs_struct)
    cpacs_main_space.set_current_variable(name_var_1, np.array([80.]))
    cpacs_main_space.set_current_variable(name_var_2, np.array([35.]))

    cpacs_subspace_1 = deepcopy(cpacs_main_space).filter(name_var_1)
    cpacs_subspace_2 = deepcopy(cpacs_main_space).filter([name_var_1, name_var_2])

    cpacs_subspace_1.set_current_variable(name_var_1, np.array([100.]))
    with pytest.raises(ValueError):
        cpacs_subspace_1.set_current_variable(name_var_2, np.array([15.]))
    cpacs_subspace_1.set_cpacs_values()

    assert cpacs_struct.get_value(name_var_1) == pytest.approx(100.)
    assert cpacs_struct.get_value(name_var_2) == pytest.approx(20.)

    cpacs_subspace_2.set_current_variable(name_var_2, np.array([18.]))
    cpacs_subspace_2.set_cpacs_values()

    assert cpacs_struct.get_value(name_var_1) == pytest.approx(80.)
    assert cpacs_struct.get_value(name_var_2) == pytest.approx(18.)

    cpacs_subspace_2.set_current_variable(name_var_1, np.array([110.]))
    cpacs_subspace_2.set_current_variable(name_var_2, np.array([17.]))
    cpacs_subspace_2.set_cpacs_values()

    assert cpacs_struct.get_value(name_var_1) == pytest.approx(110.)
    assert cpacs_struct.get_value(name_var_2) == pytest.approx(17.)
    assert cpacs_subspace_1.get_current_x([name_var_1]) == pytest.approx(100.)


def test_setting_values_with_process(cpacs_struct):
    """Test that the design space can set CPACS values with
    processed variables defined."""

    thickness_process = lambda var1, var2: 2*var1 + 3*var2**2

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                            "wing/componentSegments/componentSegment/"
                                            "structure/spars/sparSegments/"
                                            "sparSegment[1]/sparCrossSection/"
                                            "lowerCap/material/thickness",
                                            "thickness",
                                            process=thickness_process)

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/"
                                            "wings/wing/positionings/"
                                            "positioning[3]/length",
                                            "length",
                                            is_design_variable=True)

    cpacs_space = CPACSDesignSpace(cpacs_struct)
    cpacs_space.add_variable("var1", value=5.)
    cpacs_space.add_variable("var2", value=10.)
    cpacs_space.add_variable("var3", value=20.)

    cpacs_space.set_current_variable("length", np.array([55.]))

    cpacs_space.set_cpacs_values()

    assert cpacs_struct.get_value("thickness") == pytest.approx(310.)
    assert cpacs_struct.get_value("length") == pytest.approx(55.)


def test_copy_and_filter_with_processed_variables_error(cpacs_struct):
    """Test that an error is raised if arguments of a processed function
    are splitted during a design space filtering"""

    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/"
                                            "wingAirfoils/wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB",
                                            "x1", is_design_variable=True)

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                            "wing/componentSegments/componentSegment/"
                                            "structure/spars/sparSegments/"
                                            "sparSegment[1]/sparCrossSection/"
                                            "lowerCap/material/thickness",
                                            "x2",
                                            process=lambda x1, x3, x4: 3*x1+x3+2*x4)

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/"
                                            "wings/wing/positionings/"
                                            "positioning[3]/length",
                                            "x3",
                                            is_design_variable=True)

    cpacs_main_space = CPACSDesignSpace(cpacs_struct)
    cpacs_main_space.add_variable("x4", value=1.2)

    cpacs_sub_space = deepcopy(cpacs_main_space).filter(["x1", "x3"])
    cpacs_sub_space.set_current_variable("x1", np.array([10.]))
    cpacs_sub_space.set_current_variable("x3", np.array([20.]))

    with pytest.raises(ValueError):
        cpacs_sub_space.set_cpacs_values()


def test_copy_and_filter_with_processed_variables(cpacs_struct):
    """Test that processed variable are correctly set into CPACS data after
    copying and filtering variables."""

    cpacs_struct.select_variable_from_xpath("./vehicles/profiles/"
                                            "wingAirfoils/wingAirfoil[@uID='NACA_CST']/"
                                            "cst2D/upperB",
                                            "x1", is_design_variable=True)

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                            "wing/componentSegments/componentSegment/"
                                            "structure/spars/sparSegments/"
                                            "sparSegment[1]/sparCrossSection/"
                                            "lowerCap/material/thickness",
                                            "x2",
                                            process=lambda x1, x3, x4: 3*x1+x3+2*x4)

    cpacs_struct.select_variable_from_xpath("./vehicles/aircraft/model/"
                                            "wings/wing/positionings/"
                                            "positioning[3]/length",
                                            "x3",
                                            is_design_variable=True)

    cpacs_main_space = CPACSDesignSpace(cpacs_struct)
    cpacs_main_space.add_variable("x4", value=1.2)

    cpacs_sub_space = deepcopy(cpacs_main_space).filter(["x1", "x3", "x4"])
    cpacs_sub_space.set_current_variable("x1", np.array([10.]))
    cpacs_sub_space.set_current_variable("x3", np.array([20.]))
    cpacs_sub_space.set_current_variable("x4", np.array([-5.]))

    cpacs_sub_space.set_cpacs_values()

    assert cpacs_struct.get_value("x2") == pytest.approx(40.)



