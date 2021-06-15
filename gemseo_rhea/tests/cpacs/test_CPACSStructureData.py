
import pytest
import numpy as np

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSStructureData

dir_name = Path(__file__).parent.absolute()


def test_add_variable_from_empty_xpath():
    """Test a variable addition from a wrong xpath. """
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")

    with pytest.raises(ValueError, match="None element found "
         "corresponding to the current XPath."):
        cpacs.select_input_from_xpath("./titi/toto/tata", "toto")


def test_add_variable_from_several_elements_xpath():
    """Test a variable addition that returns several elements of xml tree."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")

    with pytest.raises(ValueError,
                       match="Current XPath has more than one element."):
        cpacs.select_input_from_xpath("./vehicles/aircraft/model/wings/"
        + "wing/sections/section/transformation/scaling/x",
                                      "var1")


def test_add_variable_from_xpath():
    """Test variable addition. Check that variable is found in internal dict."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")
    name_var = "x_scaling_section1_wing1"

    cpacs.select_input_from_xpath("./vehicles/aircraft/model/wings/"
    + "wing[@uID='wing1']/sections/"
    + "section[@uID='wing1section1']/transformation/scaling/x",
                                  name_var)

    assert name_var in cpacs._input_mapping._dict_elements.keys()


def test_write_xml(tmp_path):
    """Test xml writing after setting new values."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")
    name_var = "lower_naca"
    cpacs.select_input_from_xpath("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                  name_var)

    old_lower_naca_values = cpacs.get_input_values(name_var)

    new_values = np.array([4, 5, 6])
    cpacs.set_input_values(name_var, new_values)
    output_file = tmp_path / "test_output_parser.xml"
    cpacs.write_xml(output_file)

    new_cpacs = CPACSStructureData(output_file)
    new_cpacs.select_input_from_xpath("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB",
                                  name_var)
    new_lower_naca_values = new_cpacs.get_input_values("lower_naca")

    assert new_lower_naca_values == pytest.approx(new_values)
    assert old_lower_naca_values == pytest.approx(np.array([0.15, 0.1, 0.001]))


def test_add_lower_profile():
    """ Test the addition of a lower profile variable."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")
    name_var_l = "lower_naca"
    cpacs.select_lower_profile_variable(name_var_l, airfoil_ID="NACA_CST")
    assert cpacs.get_input_values(name_var_l) == pytest.approx(np.array([0.15, 0.1, 0.001]))


def test_add_upper_profile():
    """ Test the addition of an upper profile variable."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")
    name_var_u = "upper_naca"
    cpacs.select_upper_profile_variable(name_var_u, airfoil_ID="NACA_CST")
    assert cpacs.get_input_values(name_var_u) == pytest.approx(0.2)


def test_add_sweep_angle():
    """ Test the addition of a sweep angle variable."""
    cpacs = CPACSStructureData(dir_name / "test_input_wing.xml")
    name_var = "sweep"
    cpacs.select_sweep_angle_variable(name_var, wing_ID="wing1", position_ID="pos3")
    assert cpacs.get_input_values(name_var) == pytest.approx(20.)
