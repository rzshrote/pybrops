import pytest
from pybrops.breed.prot.sel.soln.MateSelectionSolution import MateSelectionSolution
from pybrops.breed.prot.sel.soln.MateSelectionSolution import check_is_MateSelectionSolution
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_isconcrete, not_raises

from .common_fixtures import *


class DummyMateSelectionSolution(MateSelectionSolution):
    def __init__(self, decn_space_xmap):
        self.decn_space_xmap = decn_space_xmap

@pytest.fixture
def selsoln(
        common_decn_space_xmap
    ):
    out = DummyMateSelectionSolution(
        decn_space_xmap = common_decn_space_xmap
    )
    yield out


################### Test class abstract/concrete properties ####################
def test_MateSelectionSolution_is_concrete():
    assert_class_isconcrete(MateSelectionSolution)

############################## Test class docstring ############################
def test_MateSelectionSolution_docstring():
    assert_class_documentation(MateSelectionSolution)

############################ Test class properties #############################
def test_MateSelectionSolution_decn_space_xmap_is_concrete():
    assert_property_isconcrete(MateSelectionSolution, "decn_space_xmap")

############################# Test class utilities #############################
def test_check_is_MateSelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_MateSelectionSolution(selsoln, "selsoln")
