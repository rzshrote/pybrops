import pytest
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import IntegerMateSelectionProtocol
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import check_is_IntegerMateSelectionProtocol
from pybrops.test.assert_python import assert_method_isabstract, assert_method_isconcrete, assert_class_documentation, assert_class_issemiabstract, not_raises

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
    assert_class_issemiabstract(IntegerMateSelectionProtocol)

############################## Test class docstring ############################
def test_IntegerMateSelectionProtocol_docstring():
    assert_class_documentation(IntegerMateSelectionProtocol)

############################## Test class methods ##############################
def test_IntegerMateSelectionProtocol_problem_is_abstract():
    assert_method_isabstract(IntegerMateSelectionProtocol, "problem")

def test_IntegerMateSelectionProtocol_sosolve_is_concrete():
    assert_method_isconcrete(IntegerMateSelectionProtocol, "sosolve")

def test_IntegerMateSelectionProtocol_mosolve_is_concrete():
    assert_method_isconcrete(IntegerMateSelectionProtocol, "mosolve")

def test_IntegerMateSelectionProtocol_select_is_concrete():
    assert_method_isconcrete(IntegerMateSelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_IntegerMateSelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_IntegerMateSelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_IntegerMateSelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_IntegerMateSelectionProtocol(None, "selsoln")
