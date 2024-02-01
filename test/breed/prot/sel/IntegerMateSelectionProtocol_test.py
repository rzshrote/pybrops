import pytest
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import IntegerMateSelectionProtocol
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import check_is_IntegerMateSelectionProtocol
from pybrops.test.assert_python import assert_abstract_method, assert_concrete_method, assert_docstring, assert_semiabstract_class, not_raises

class DummyIntegerMateSelectionProtocol(IntegerMateSelectionProtocol):
    def __init__(self):
        pass
    def problem(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        pass

@pytest.fixture
def selsoln():
    out = DummyIntegerMateSelectionProtocol()
    yield out

################### Test class abstract/concrete properties ####################
def test_IntegerMateSelectionProtocol_is_semiabstract():
    assert_semiabstract_class(IntegerMateSelectionProtocol)

############################## Test class docstring ############################
def test_IntegerMateSelectionProtocol_docstring():
    assert_docstring(IntegerMateSelectionProtocol)

############################## Test class methods ##############################
def test_IntegerMateSelectionProtocol_problem_is_abstract():
    assert_abstract_method(IntegerMateSelectionProtocol, "problem")

def test_IntegerMateSelectionProtocol_sosolve_is_concrete():
    assert_concrete_method(IntegerMateSelectionProtocol, "sosolve")

def test_IntegerMateSelectionProtocol_mosolve_is_concrete():
    assert_concrete_method(IntegerMateSelectionProtocol, "mosolve")

def test_IntegerMateSelectionProtocol_select_is_concrete():
    assert_concrete_method(IntegerMateSelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_IntegerMateSelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_IntegerMateSelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_IntegerMateSelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_IntegerMateSelectionProtocol(None, "selsoln")
