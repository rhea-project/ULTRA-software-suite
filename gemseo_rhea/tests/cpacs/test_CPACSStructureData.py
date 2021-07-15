
import pytest
import numpy as np

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSStructureData

dir_name = Path(__file__).parent.absolute()


@pytest.fixture
def cpacs():
    return CPACSStructureData(dir_name / "test_input_wing.xml")


def test_add_variable_from_empty_xpath(cpacs):
    """Test a variable addition from a wrong xpath. """

    with pytest.raises(ValueError, match="None element found "
         "corresponding to the current XPath."):
        cpacs.select_variable_from_xpath("./titi/toto/tata", "toto")


def test_add_variable_from_several_elements_xpath(cpacs):
    """Test a variable addition that returns several elements of xml tree."""

    with pytest.raises(ValueError,
                       match="Current XPath has more than one element."):
        cpacs.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                         + "wing/sections/section/transformation/scaling/x",
                                      "var1")


def test_add_variable_from_xpath(cpacs):
    """Test variable addition. Check that variable is found in internal dict."""
    name_var = "x_scaling_section1_wing1"

    cpacs.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                     + "wing[@uID='wing1']/sections/"
                                     + "section[@uID='wing1section1']/"
                                       "transformation/scaling/x",
                                     name_var)

    assert cpacs.get_value(name_var) == pytest.approx(1.5)


def test_write_xml(tmp_path, cpacs):
    """Test xml writing after setting new values."""
    name_var = "lower_naca"
    cpacs.select_variable_from_xpath("./vehicles/profiles/"
                                     + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                     name_var)

    old_lower_naca_values = cpacs.get_value(name_var)

    new_values = np.array([4, 5, 6])
    cpacs.set_value(name_var, value=new_values)
    output_file = tmp_path / "test_output_parser.xml"
    cpacs.write_xml(output_file)

    new_cpacs = CPACSStructureData(output_file)
    new_cpacs.select_variable_from_xpath("./vehicles/profiles/"
                                         + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                         name_var)
    new_lower_naca_values = new_cpacs.get_value("lower_naca")

    assert new_lower_naca_values == pytest.approx(new_values)
    assert old_lower_naca_values == pytest.approx(np.array([0.15, 0.1, 0.001]))


def test_sub_mapping(cpacs):
    """Test that sub mapping is correctly created."""
    cpacs.select_lower_profile_variable("x1", airfoil_ID="NACA_CST")
    cpacs.select_upper_profile_variable("x2", airfoil_ID="NACA_CST")
    cpacs.select_sweep_angle_variable("x3", wing_ID="wing1", position_ID="pos3")
    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                      "componentSegment/structure/"
                                      "ribsDefinitions/ribsDefinition/"
                                      "ribsPositioning/ribRotation/z",
                                      "x4")

    sub_mapping = cpacs.get_sub_mapping(["x2", "x4"])

    assert sub_mapping["x2"] == pytest.approx(0.2)
    assert sub_mapping["x4"] == pytest.approx(90.)


def test_iterator_with_selected_var(cpacs):
    """Test iterating over the cpacs structure."""
    cpacs.select_lower_profile_variable("x1", airfoil_ID="NACA_CST",
                                        is_design_variable=True)

    cpacs.select_upper_profile_variable("x2", airfoil_ID="NACA_CST",
                                        is_design_variable=True,
                                        lower_bound=-10, upper_bound=10)

    cpacs.select_sweep_angle_variable("x3", wing_ID="wing1", position_ID="pos3")
    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                      "componentSegment/structure/"
                                      "ribsDefinitions/ribsDefinition/"
                                      "ribsPositioning/ribRotation/z",
                                      "x4", is_design_variable=True,
                                     lower_bound=-5.2, upper_bound=9.6)

    variables = {}
    variables_bounds = {}
    for var, values in cpacs:
        variables.update({var: values})
        variables_bounds.update({var: cpacs.get_bounds(var)})

    assert len(variables) == 3
    assert variables_bounds["x1"] == (None, None)
    assert variables["x1"] == pytest.approx(np.array([0.15, 0.1, 0.001]))
    assert variables_bounds["x2"] == pytest.approx((-10, 10))
    assert variables["x2"] == pytest.approx(0.2)
    assert variables_bounds["x4"] == pytest.approx((-5.2, 9.6))
    assert variables["x4"] == pytest.approx(90)


def test_iterator_with_no_var(cpacs):
    """Test iterating over the cpacs structure."""
    cpacs.select_lower_profile_variable("x1", airfoil_ID="NACA_CST")
    cpacs.select_upper_profile_variable("x2", airfoil_ID="NACA_CST")

    variables = []
    for var, values in cpacs:
        variables.append(var)

    assert len(variables) == 0


def test_processed_variable_setting_values(cpacs):
    """Test that setting values to processed variable works fine."""
    rib_rot_update = lambda v1, v2, v3: 2*v1 + v2**2 - 3*v3

    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                      "componentSegment/structure/"
                                      "ribsDefinitions/ribsDefinition/"
                                      "ribsPositioning/ribRotation/z",
                                      "x4", process=rib_rot_update)

    cpacs.set_value("x4", v1=2, v2=1, v3=4)

    assert cpacs.get_value("x4") == pytest.approx(-7)


def test_processed_variable_arguments(cpacs):
    """Test that process arguments are returned."""
    rib_rot_update = lambda v1, v2, v3: 2*v1 + v2**2 - 3*v3

    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                      "componentSegment/structure/"
                                      "ribsDefinitions/ribsDefinition/"
                                      "ribsPositioning/ribRotation/z",
                                      "x4", process=rib_rot_update)

    cpacs.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                     "wing[@uID='wing1']/sections/"
                                     "section[@uID='wing1section1']/"
                                     "transformation/scaling/x",
                                     "x3", process=lambda t: 2*t)

    assert list(cpacs.get_process_arguments("x4")) == ["v1", "v2", "v3"]


def test_processed_variable_error_with_design_space(cpacs):
    """Test that a processed variable cannot be a design variable."""
    rib_rot_update = lambda v1, v2, v3: 2*v1 + v2**2 - 3*v3

    with pytest.raises(ValueError, match="CPACS variable cannot be a design variable"
                                 "and a processed variable in the same time."):
        cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                         "model/wings/wing/componentSegments/"
                                          "componentSegment/structure/"
                                          "ribsDefinitions/ribsDefinition/"
                                          "ribsPositioning/ribRotation/z",
                                          "x4", process=rib_rot_update,
                                         is_design_variable=True)


def test_design_variables_property(cpacs):
    """Test that design variables are returned."""

    cpacs.select_lower_profile_variable("x1", airfoil_ID="NACA_CST",
                                        is_design_variable=True)

    cpacs.select_upper_profile_variable("x2", airfoil_ID="NACA_CST",
                                        is_design_variable=True)

    cpacs.select_sweep_angle_variable("x3", wing_ID="wing1", position_ID="pos3",
                                      is_design_variable=True)

    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                     "componentSegment/structure/"
                                     "ribsDefinitions/ribsDefinition/"
                                     "ribsPositioning/ribRotation/z",
                                     "x4",
                                     process=lambda x1:2*x1)

    assert len(cpacs.design_variables) == 3
    assert "x1" in cpacs.design_variables
    assert "x2" in cpacs.design_variables
    assert "x3" in cpacs.design_variables


def test_processed_variables_property(cpacs):
    """Test that processed variables are returned."""

    cpacs.select_variable_from_xpath("./vehicles/profiles/wingAirfoils/"
                                     "wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                     "x1",
                                     process=lambda x3:3*x3)

    cpacs.select_upper_profile_variable("x2", airfoil_ID="NACA_CST",
                                        is_design_variable=True)

    cpacs.select_sweep_angle_variable("x3", wing_ID="wing1", position_ID="pos3",
                                      is_design_variable=True)

    cpacs.select_variable_from_xpath("./vehicles/aircraft/"
                                     "model/wings/wing/componentSegments/"
                                     "componentSegment/structure/"
                                     "ribsDefinitions/ribsDefinition/"
                                     "ribsPositioning/ribRotation/z",
                                     "x4",
                                     process=lambda x2:2*x2)

    assert len(cpacs.processed_variables) == 2
    assert "x1" in cpacs.processed_variables
    assert "x4" in cpacs.processed_variables


def test_add_lower_profile(cpacs):
    """ Test the addition of a lower profile variable."""
    name_var_l = "lower_naca"
    cpacs.select_lower_profile_variable(name_var_l, airfoil_ID="NACA_CST")
    assert cpacs.get_value(name_var_l) == pytest.approx(np.array([0.15, 0.1, 0.001]))


def test_add_upper_profile(cpacs):
    """ Test the addition of an upper profile variable."""
    name_var_u = "upper_naca"
    cpacs.select_upper_profile_variable(name_var_u, airfoil_ID="NACA_CST")
    assert cpacs.get_value(name_var_u) == pytest.approx(0.2)


def test_add_sweep_angle(cpacs):
    """ Test the addition of a sweep angle variable."""
    name_var = "sweep"
    cpacs.select_sweep_angle_variable(name_var, wing_ID="wing1", position_ID="pos3")
    assert cpacs.get_value(name_var) == pytest.approx(20.)


def test_add_twist_angle(cpacs):
    """Test the selection of a twist angle variable."""
    name_var = "twist"
    cpacs.select_twist_angle(name_var, "wing1", "wing1section2")
    assert cpacs.get_value(name_var) == pytest.approx(1.2)