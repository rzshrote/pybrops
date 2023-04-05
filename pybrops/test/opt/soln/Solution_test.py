import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.Solution import Solution, check_is_Solution

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield Solution()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(Solution)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(Solution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ndecn_is_abstract():
    assert_abstract_property(Solution, "ndecn")

def test_decn_space_is_abstract():
    assert_abstract_property(Solution, "decn_space")

def test_decn_space_lower_is_abstract():
    assert_abstract_property(Solution, "decn_space_lower")

def test_decn_space_upper_is_abstract():
    assert_abstract_property(Solution, "decn_space_upper")

def test_nobj_is_abstract():
    assert_abstract_property(Solution, "nobj")

def test_obj_wt_is_abstract():
    assert_abstract_property(Solution, "obj_wt")

def test_nineqcv_is_abstract():
    assert_abstract_property(Solution, "nineqcv")

def test_ineqcv_wt_is_abstract():
    assert_abstract_property(Solution, "ineqcv_wt")

def test_neqcv_is_abstract():
    assert_abstract_property(Solution, "neqcv")

def test_eqcv_wt_is_abstract():
    assert_abstract_property(Solution, "eqcv_wt")

def test_soln_decn_is_abstract():
    assert_abstract_property(Solution, "soln_decn")

def test_soln_obj_is_abstract():
    assert_abstract_property(Solution, "soln_obj")

def test_soln_ineqcv_is_abstract():
    assert_abstract_property(Solution, "soln_ineqcv")

def test_soln_eqcv_is_abstract():
    assert_abstract_property(Solution, "soln_eqcv")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_Solution_is_concrete():
    assert_concrete_function(check_is_Solution)

def test_check_is_Solution(prob):
    with not_raises(TypeError):
        check_is_Solution(prob, "prob")
    with pytest.raises(TypeError):
        check_is_Solution(None, "prob")
