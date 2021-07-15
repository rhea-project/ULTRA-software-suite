
import pytest
import numpy as np

from gemseo_rhea.elementary_bilevel.cpacs import CPACS_MR_SBW
from gemseo_rhea.elementary_bilevel.design_space import DesignSpace_MR_SBW


@pytest.fixture
def design_space():
    cpacs = CPACS_MR_SBW()
    cpacs.set_planform_design_variables()
    cpacs.set_planform_processed_variables()
    space = DesignSpace_MR_SBW(cpacs)
    return space


def test_planform_variables_definition(design_space):
    """Test that planform variables are correctly set."""
    design_space.add_planform_variables()

    assert len(design_space.variables_names) == 6

    assert "span" in design_space.variables_names
    assert "taper_ratio" in design_space.variables_names
    assert "root_chord" in design_space.variables_names
    assert "twist_kink" in design_space.variables_names
    assert "twist_tip" in design_space.variables_names
    assert "sweep_section_1" in design_space.variables_names


def test_planform_variables_setting_values(design_space):
    """Test that setting values to design variables update
    processed CPACS variables"""

    design_space.add_planform_variables()

    design_space.set_current_variable("span", np.array([20]))
    design_space.set_cpacs_values()

    assert design_space._cpacs_structure.get_value("section_length_1") == \
           pytest.approx(10.)
    assert design_space._cpacs_structure.get_value("section_length_2") == \
           pytest.approx(10.)


