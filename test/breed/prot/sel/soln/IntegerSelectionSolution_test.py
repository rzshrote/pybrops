import pytest
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.breed.prot.sel.soln.IntegerSelectionSolution import IntegerSelectionSolution
from pybrops.breed.prot.sel.soln.IntegerSelectionSolution import check_is_IntegerSelectionSolution
from pybrops.test.assert_python import assert_class_documentation, assert_class_isconcrete, not_raises


class DummyIntegerSelectionSolution(IntegerSelectionSolution):
    def __init__(self):
        pass
    

@pytest.fixture
def selsoln():
    out = DummyIntegerSelectionSolution()
    yield out


################### Test class abstract/concrete properties ####################
def test_IntegerSelectionSolution_is_concrete():
    assert_class_isconcrete(IntegerSelectionSolution)

############################## Test class docstring ############################
def test_IntegerSelectionSolution_docstring():
    assert_class_documentation(IntegerSelectionSolution)

############################# Test class utilities #############################
def test_check_is_IntegerSelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_IntegerSelectionSolution(selsoln, "selsoln")
