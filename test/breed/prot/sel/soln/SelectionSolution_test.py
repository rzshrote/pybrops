import pytest
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.breed.prot.sel.soln.SelectionSolution import check_is_SelectionSolution
from pybrops.test.assert_python import assert_class_documentation, assert_class_isconcrete, not_raises


class DummySelectionSolution(SelectionSolution):
    def __init__(self):
        pass
    

@pytest.fixture
def selsoln():
    out = DummySelectionSolution()
    yield out


################### Test class abstract/concrete properties ####################
def test_SelectionSolution_is_concrete():
    assert_class_isconcrete(SelectionSolution)

############################## Test class docstring ############################
def test_SelectionSolution_docstring():
    assert_class_documentation(SelectionSolution)

############################# Test class utilities #############################
def test_check_is_SelectionSolution(selsoln):
    with not_raises(Exception):
        check_is_SelectionSolution(selsoln, "selsoln")
