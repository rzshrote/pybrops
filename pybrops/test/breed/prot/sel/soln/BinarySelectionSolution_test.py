import pytest
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.breed.prot.sel.soln.BinarySelectionSolution import BinarySelectionSolution
from pybrops.breed.prot.sel.soln.BinarySelectionSolution import check_is_BinarySelectionSolution
from pybrops.test.assert_python import assert_docstring, assert_concrete_class, not_raises


class DummyBinarySelectionSolution(BinarySelectionSolution):
    def __init__(self):
        pass
    

@pytest.fixture
def selsoln():
    out = DummyBinarySelectionSolution()
    yield out


################### Test class abstract/concrete properties ####################
def test_BinarySelectionSolution_is_concrete():
    assert_concrete_class(BinarySelectionSolution)

############################## Test class docstring ############################
def test_BinarySelectionSolution_docstring():
    assert_docstring(BinarySelectionSolution)

############################# Test class utilities #############################
def test_check_is_BinarySelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_BinarySelectionSolution(selsoln, "selsoln")
