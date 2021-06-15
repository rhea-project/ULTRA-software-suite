
import pytest

from gemseo_rhea.cpacs import CPACSVariable


def test_xpath_template_1():
    var = CPACSVariable("./toto[@uID='{airfoil_ID}']/titi",
        ["airfoil_ID"])
    assert var.xpath_template == "./toto[@uID='{airfoil_ID}']/titi"


def test_list_id_template():
    var = CPACSVariable("./toto[@uID='{airfoil_ID}']/titi",
        ["airfoil_ID"])
    assert var.list_id_template == ["airfoil_ID"]


def test_get_xpath():
    var = CPACSVariable("./toto[@uID='{airfoil_ID}']/titi",
        ["airfoil_ID"])
    xpath = var.get_xpath(["NACA_CST"])
    assert xpath == "./toto[@uID='NACA_CST']/titi"


def test_get_xpath_3_vars():
    var = CPACSVariable("./toto[@uID='{v1}']/titi"
                        + "/tutu[@uID='{v2}']/tete[@uID='{v3}'",
        ["v1", "v2", "v3"])
    xpath = var.get_xpath(["WING_1", "SEC_1", "RIB_3"])
    assert xpath == "./toto[@uID='WING_1']/titi/tutu[@uID='SEC_1']/tete[@uID='RIB_3'"


def test_get_xpath_error():
    var = CPACSVariable("./toto[@uID='{v1}']/titi"
                        + "/tutu[@uID='{v2}']/tete[@uID='{v3}'",
        ["v1", "v2", "v3"])
    with pytest.raises(ValueError,
        match="The number of values in list_id_values should be 3 while it is 2."):
        var.get_xpath(["WING_1", "SEC_1"])


