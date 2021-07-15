
import pytest
from pathlib import Path

from .RazorParser import RazorParser


@pytest.fixture
def razor_parser():
    name_file = Path(__file__).parent / "tutorial_naca0012_vol.cfg"
    return RazorParser(name_file)


def test_variables(razor_parser):
    """Check the name of variables."""
    assert razor_parser.variables == ["TIME"]


def test_fields(razor_parser):
    """Check the name of fields."""
    assert razor_parser.fields == ["X-Momentum", "Y-Momentum"]

def test_database(razor_parser):
    assert razor_parser._RazorParser__database == "naca0012_test_volume.h5"
    

def test_build_target(razor_parser):
    """Check that target list is correctly written."""
    target = razor_parser._RazorParser__build_target_list({"TIME":0.5}, "output.csv")
    assert target == "TARGET_LIST=0.5,output.csv\n"


def test_target_file_generation(tmp_path, razor_parser):
    """Check that target file is built."""
    new_input = tmp_path / "new_input_test.csv"
    razor_parser.generate_target_file({"TIME": 0.3}, output_target_file="out_test.csv",
                                      generated_file=new_input)

    f=open(new_input, "r")
    line = f.readline()
    while "TARGET_LIST=" not in line:
        line = f.readline()
    assert line.strip() == "TARGET_LIST=0.3,out_test.csv"


def test_output_parsing_SU2_CSV(razor_parser):
    """Check that the ouputs is correctly parsed with SU2 format."""
    name_file = Path(__file__).parent / "razor_output_for_test.csv"
    result = razor_parser.parse_output(name_file)
    assert len(result) == 2
    assert "X-Momentum" in result.keys()
    assert "Y-Momentum" in result.keys()
