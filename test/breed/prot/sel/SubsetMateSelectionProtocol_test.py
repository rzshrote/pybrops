import pytest
from pybrops.breed.prot.sel.SubsetMateSelectionProtocol import SubsetMateSelectionProtocol
from pybrops.breed.prot.sel.SubsetMateSelectionProtocol import check_is_SubsetMateSelectionProtocol
from pybrops.test.assert_python import assert_abstract_method, assert_concrete_method, assert_docstring, assert_semiabstract_class, not_raises

class DummySubsetMateSelectionProtocol(SubsetMateSelectionProtocol):
    def __init__(self):
        pass
    def problem(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        pass

@pytest.fixture
def selsoln():
    out = DummySubsetMateSelectionProtocol()
    yield out

################### Test class abstract/concrete properties ####################
def test_SubsetMateSelectionProtocol_is_semiabstract():
    assert_semiabstract_class(SubsetMateSelectionProtocol)

############################## Test class docstring ############################
def test_SubsetMateSelectionProtocol_docstring():
    assert_docstring(SubsetMateSelectionProtocol)

############################## Test class methods ##############################
def test_SubsetMateSelectionProtocol_problem_is_abstract():
    assert_abstract_method(SubsetMateSelectionProtocol, "problem")

def test_SubsetMateSelectionProtocol_sosolve_is_concrete():
    assert_concrete_method(SubsetMateSelectionProtocol, "sosolve")

def test_SubsetMateSelectionProtocol_mosolve_is_concrete():
    assert_concrete_method(SubsetMateSelectionProtocol, "mosolve")

def test_SubsetMateSelectionProtocol_select_is_concrete():
    assert_concrete_method(SubsetMateSelectionProtocol, "select")

############################# Test class utilities #############################
def test_check_is_SubsetMateSelectionProtocol(selsoln):
    with not_raises(Exception):
        check_is_SubsetMateSelectionProtocol(selsoln, "selsoln")
    with pytest.raises(TypeError):
        check_is_SubsetMateSelectionProtocol(object(), "selsoln")
    with pytest.raises(TypeError):
        check_is_SubsetMateSelectionProtocol(None, "selsoln")
