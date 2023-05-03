import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.SolutionType import SolutionType, check_is_SolutionType

################################################################################
################################ Test fixtures #################################
################################################################################

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(SolutionType)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init():
    with pytest.raises(TypeError):
        soln = SolutionType()

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ndecn_is_abstract():
    assert_abstract_property(SolutionType, "ndecn")

def test_decn_space_is_abstract():
    assert_abstract_property(SolutionType, "decn_space")

def test_decn_space_lower_is_abstract():
    assert_abstract_property(SolutionType, "decn_space_lower")

def test_decn_space_upper_is_abstract():
    assert_abstract_property(SolutionType, "decn_space_upper")

def test_nobj_is_abstract():
    assert_abstract_property(SolutionType, "nobj")

def test_obj_wt_is_abstract():
    assert_abstract_property(SolutionType, "obj_wt")

def test_nineqcv_is_abstract():
    assert_abstract_property(SolutionType, "nineqcv")

def test_ineqcv_wt_is_abstract():
    assert_abstract_property(SolutionType, "ineqcv_wt")

def test_neqcv_is_abstract():
    assert_abstract_property(SolutionType, "neqcv")

def test_eqcv_wt_is_abstract():
    assert_abstract_property(SolutionType, "eqcv_wt")

def test_soln_decn_is_abstract():
    assert_abstract_property(SolutionType, "soln_decn")

def test_soln_obj_is_abstract():
    assert_abstract_property(SolutionType, "soln_obj")

def test_soln_ineqcv_is_abstract():
    assert_abstract_property(SolutionType, "soln_ineqcv")

def test_soln_eqcv_is_abstract():
    assert_abstract_property(SolutionType, "soln_eqcv")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SolutionType_is_concrete():
    assert_concrete_function(check_is_SolutionType)

def test_check_is_SolutionType():
    with pytest.raises(TypeError):
        check_is_SolutionType(None, "soln")
