import pytest
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.breed.prot.sel.soln.RealSelectionSolution import RealSelectionSolution
from pybrops.breed.prot.sel.soln.RealSelectionSolution import check_is_RealSelectionSolution
from pybrops.test.assert_python import assert_class_documentation, assert_class_isconcrete, not_raises


class DummyRealSelectionSolution(RealSelectionSolution):
    def __init__(self):
        pass
    

@pytest.fixture
def selsoln():
    out = DummyRealSelectionSolution()
    yield out


################### Test class abstract/concrete properties ####################
def test_RealSelectionSolution_is_concrete():
    assert_class_isconcrete(RealSelectionSolution)

############################## Test class docstring ############################
def test_RealSelectionSolution_docstring():
    assert_class_documentation(RealSelectionSolution)

############################# Test class utilities #############################
def test_check_is_RealSelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_RealSelectionSolution(selsoln, "selsoln")
