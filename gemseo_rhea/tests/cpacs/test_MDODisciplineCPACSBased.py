
import numpy as np
import pytest

from gemseo.utils.py23_compat import Path
from gemseo_rhea.cpacs import MDODisciplineCPACSBased
from gemseo_rhea.cpacs import CPACSStructureData


class TestDiscipline(MDODisciplineCPACSBased):

    def _run(self):
        x1 = self.get_inputs_by_name("x1")
        x2 = self.get_inputs_by_name("x2")
        y1 = x1**2 + 2*x2**3
        self.local_data["y1"] = y1

    def _compute_jacobian(self, inputs=None, outputs=None):
        self._init_jacobian(with_zeros=True)
        x1 = self.get_inputs_by_name("x1")
        x2 = self.get_inputs_by_name("x2")
        self.jac["y1"] = {}
        self.jac["y1"]["x1"] = np.atleast_2d(2*x1)
        self.jac["y1"]["x2"] = np.atleast_2d(6*x2**2)


@pytest.fixture
def discipline():
    input_cpacs_file = Path(__file__).parent / "test_input_wing.xml"
    cpacs = CPACSStructureData(input_cpacs_file)

    cpacs.select_variable_from_xpath(xpath="./vehicles/aircraft/model/"
                                                "wings/wing[@uID='wing1']/"
                                                "positionings/positioning[@uID='pos3']/"
                                                "sweepAngle",
                                          name="x1",
                                          is_design_variable=True,
                                          lower_bound=-20.,
                                          upper_bound=20.)

    cpacs.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                     "wing[@uID='wing1']/sections/"
                                     "section[@uID='wing1section1']/"
                                     "transformation/scaling/x",
                                     name="x2",
                                     is_design_variable=True,
                                     lower_bound=-5.,
                                     upper_bound=5.)

    cpacs.select_variable_from_xpath("./vehicles/aircraft/model/wings/"
                                          "wing/componentSegments/componentSegment/"
                                          "structure/spars/sparPositions/"
                                          "sparPosition[3]/sparPositionEtaXsi/"
                                          "xsi", "y1")

    in_mapping = cpacs.get_sub_mapping(["x1", "x2"])
    out_mapping = cpacs.get_sub_mapping(["y1"])

    discipline = TestDiscipline(in_mapping, out_mapping)

    return discipline


def test_discipline_building(discipline):
    """Test that discipline grammar is correctly set."""
    assert discipline.input_grammar.is_data_name_existing("x1")
    assert discipline.input_grammar.is_data_name_existing("x2")
    assert discipline.output_grammar.is_data_name_existing("y1")


def test_discipline_exec_default_values(discipline):
    """Test that discipline is correctly executed with default values."""
    result = discipline.execute()
    assert result["x1"] == pytest.approx(20.)
    assert result["x2"] == pytest.approx(1.5)
    assert result["y1"] == pytest.approx(406.75)


def test_setting_cpacs_values(discipline):
    """Test that cpacs values are correctly set be the discipline."""
    result = discipline.execute({"x1": np.array([2.]), "x2": np.array([1.])})
    discipline.set_cpacs_values()
    assert result["y1"] == pytest.approx(6.)
    assert discipline._MDODisciplineCPACSBased__input_cpacs_mapping["x1"] == \
           pytest.approx(2.)
    assert discipline._MDODisciplineCPACSBased__input_cpacs_mapping["x2"] == \
           pytest.approx(1.)
    assert discipline._MDODisciplineCPACSBased__output_cpacs_mapping["y1"] == \
           pytest.approx(6.)

