import pytest
from pybrops.breed.prot.sel.RealMateSelectionProtocol import RealMateSelectionProtocol
from pybrops.breed.prot.sel.RealMateSelectionProtocol import check_is_RealMateSelectionProtocol
from pybrops.test.assert_python import assert_method_isabstract, assert_method_isconcrete, assert_class_documentation, assert_class_issemiabstract, not_raises

class DummyRealMateSelectionProtocol(RealMateSelectionProtocol):
    def __init__(self):
        pass
    def problem(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        pass

@pytest.fixture
def selsoln():
    out = DummyRealMateSelectionProtocol()
    yield out

################### Test class abstract/concrete properties ####################
def test_RealMateSelectionProtocol_is_semiabstract():
    assert_class_issemiabstract(RealMateSelectionProtocol)

############################## Test class docstring ############################
def test_RealMateSelectionProtocol_docstring():
    assert_class_documentation(RealMateSelectionProtocol)

############################## Test class methods ##############################
def test_RealMateSelectionProtocol_problem_is_abstract():
    assert_method_isabstract(RealMateSelectionProtocol, "problem")

def test_RealMateSelectionProtocol_sosolve_is_concrete():
    assert_method_isconcrete(RealMateSelectionProtocol, "sosolve")

def test_RealMateSelectionProtocol_mosolve_is_concrete():
    assert_method_isconcrete(RealMateSelectionProtocol, "mosolve")

def test_RealMateSelectionProtocol_select_is_concrete():
    assert_method_isconcrete(RealMateSelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_RealMateSelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_RealMateSelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_RealMateSelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_RealMateSelectionProtocol(None, "selsoln")
