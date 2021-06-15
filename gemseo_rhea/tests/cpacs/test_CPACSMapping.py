
import pytest
import numpy as np
import xml.etree.ElementTree as ET

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import CPACSMapping

dir_name = Path(__file__).parent.absolute()


def get_xml_element(xpath):
    """Function that get an element in xml tree from xpath"""
    tree = ET.parse(dir_name / "test_input_wing.xml")
    elements = tree.findall(xpath)
    return elements[0]


def test_get_single_value_element():
    """"Test that an added variable has the right value."""
    elm = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB")

    cpacs = CPACSMapping()
    name_var = "upper_naca"
    cpacs.add_xml_element(name_var, elm)

    assert cpacs[name_var] == pytest.approx(0.2)


def test_get_multiple_values_element():
    """Test that the added variable has the right values in case
    of multiple size."""
    elm = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB")

    cpacs = CPACSMapping()
    name_var = "lower_naca"
    cpacs.add_xml_element(name_var, elm)

    assert (cpacs[name_var] ==
            pytest.approx(np.array([0.15, 0.1, 0.001])))


def test_get_non_existing_var():
    """Test that an error is raised when getting a non existing var"""
    cpacs = CPACSMapping()
    with pytest.raises(KeyError):
        cpacs["my_var"]


def test_set_values_with_wrong_size():
    """Test values setting with the wrong size."""
    elm = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB")

    cpacs = CPACSMapping()
    name_var = "lower_naca"
    cpacs.add_xml_element(name_var, elm)

    with pytest.raises(ValueError, match="Variable {}. Prescribed"
                       " values has the wrong size.".format(name_var)):
        cpacs[name_var] = np.array([10., 20.])


def test_set_values_with_good_size_3():
    """Test values setting of size 3."""
    elm = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB")

    cpacs = CPACSMapping()
    name_var = "lower_naca"
    cpacs.add_xml_element(name_var, elm)

    new_values = np.array([11, 22, 33])
    cpacs[name_var] = new_values

    new_text = cpacs._dict_elements[name_var].text

    assert new_text == "11;22;33"


def test_set_values_with_good_size_1():
    """Test value setting of size 1."""
    elm = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB")

    cpacs = CPACSMapping()
    name_var = "upper_naca"
    cpacs.add_xml_element(name_var, elm)

    new_values = np.array([999])
    cpacs[name_var] = new_values

    new_text = cpacs._dict_elements[name_var].text

    assert new_text == "999"


def test_print():
    """Test xml printing after loading variables."""
    elm_up = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB")
    elm_low = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB")

    cpacs = CPACSMapping()
    name_var_l = "lower_naca"
    cpacs.add_xml_element(name_var_l, elm_low)
    name_var_u = "upper_naca"
    cpacs.add_xml_element(name_var_u, elm_up)

    ref = "- " + name_var_l + "\n" + "\t0.15;0.1;0.001\n"
    ref+= "- " + name_var_u + "\n" + "\t0.2\n"

    assert cpacs.__str__() == ref


def test_get_dict_with_values():
    elm_up = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/upperB")
    elm_low = get_xml_element("./vehicles/profiles/"
    + "wingAirfoils/wingAirfoil[@uID='NACA_CST']/cst2D/lowerB")

    cpacs = CPACSMapping()
    name_var_l = "lower_naca"
    cpacs.add_xml_element(name_var_l, elm_low)
    name_var_u = "upper_naca"
    cpacs.add_xml_element(name_var_u, elm_up)

    res_d = cpacs.get_dict_with_values()

    assert res_d["lower_naca"] == pytest.approx(np.array([0.15, 0.1, 0.001]))
    assert res_d["upper_naca"] == pytest.approx(np.array([0.2]))
