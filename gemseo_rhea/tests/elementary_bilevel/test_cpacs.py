import pytest

from gemseo_rhea.elementary_bilevel.cpacs import CPACS_MR_SBW


@pytest.fixture
def sbw_cpacs():
    return CPACS_MR_SBW()


def test_planform_variables(sbw_cpacs):
    """Test planform values selected in CPACS"""
    sbw_cpacs.set_planform_design_variables()

    assert sbw_cpacs.get_value("twist_kink") == pytest.approx(-1.5)
    assert sbw_cpacs.get_value("twist_tip") == pytest.approx(-3.5)
    assert sbw_cpacs.get_value("root_chord") == pytest.approx(2.9848)
    assert sbw_cpacs.get_value("sweep_section_1") == pytest.approx(12.5)


def test_planform_processed_variables(sbw_cpacs):
    """Test planform variables that must be processed."""
    sbw_cpacs.set_planform_processed_variables()

    sbw_cpacs.set_value("tip_chord", root_chord=10., taper_ratio=0.5)
    assert sbw_cpacs.get_value("tip_chord") == pytest.approx(5.)

    sbw_cpacs.set_value("section_length_1", span=30.)
    assert sbw_cpacs.get_value("section_length_1") == pytest.approx(15.)

    sbw_cpacs.set_value("section_length_2", span=32.)
    assert sbw_cpacs.get_value("section_length_2") == pytest.approx(16.)

    sbw_cpacs.set_value("sweep_section_2", sweep_section_1=13.5)
    assert sbw_cpacs.get_value("sweep_section_2") == pytest.approx(13.5)