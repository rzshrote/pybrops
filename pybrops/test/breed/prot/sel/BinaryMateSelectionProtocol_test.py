import pytest
from pybrops.breed.prot.sel.BinaryMateSelectionProtocol import BinaryMateSelectionProtocol
from pybrops.breed.prot.sel.BinaryMateSelectionProtocol import check_is_BinaryMateSelectionProtocol
from pybrops.test.assert_python import assert_abstract_method, assert_concrete_method, assert_docstring, assert_semiabstract_class, not_raises

class DummyBinaryMateSelectionProtocol(BinaryMateSelectionProtocol):
    def __init__(self):
        pass
    def problem(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        pass

@pytest.fixture
def selsoln():
    out = DummyBinaryMateSelectionProtocol()
    yield out

################### Test class abstract/concrete properties ####################
def test_BinaryMateSelectionProtocol_is_semiabstract():
    assert_semiabstract_class(BinaryMateSelectionProtocol)

############################## Test class docstring ############################
def test_BinaryMateSelectionProtocol_docstring():
    assert_docstring(BinaryMateSelectionProtocol)

############################## Test class methods ##############################
def test_BinaryMateSelectionProtocol_problem_is_abstract():
    assert_abstract_method(BinaryMateSelectionProtocol, "problem")

def test_BinaryMateSelectionProtocol_sosolve_is_concrete():
    assert_concrete_method(BinaryMateSelectionProtocol, "sosolve")

def test_BinaryMateSelectionProtocol_mosolve_is_concrete():
    assert_concrete_method(BinaryMateSelectionProtocol, "mosolve")

def test_BinaryMateSelectionProtocol_select_is_concrete():
    assert_concrete_method(BinaryMateSelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_BinaryMateSelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_BinaryMateSelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_BinaryMateSelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_BinaryMateSelectionProtocol(None, "selsoln")
