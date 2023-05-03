import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.prob.ProblemType import ProblemType, check_is_ProblemType

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield ProblemType()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(ProblemType)

################################################################################
########################### Test concrete properties ###########################
################################################################################
def test_n_var_is_concrete():
    assert_concrete_property(ProblemType, "n_var")

def test_n_obj_is_concrete():
    assert_concrete_property(ProblemType, "n_obj")

def test_n_ieq_constr_is_concrete():
    assert_concrete_property(ProblemType, "n_ieq_constr")

def test_n_eq_constr_is_concrete():
    assert_concrete_property(ProblemType, "n_eq_constr")

def test_xl_is_concrete():
    assert_concrete_property(ProblemType, "xl")

def test_xu_is_concrete():
    assert_concrete_property(ProblemType, "xu")

def test_vtype_is_concrete():
    assert_concrete_property(ProblemType, "vtype")

def test_vars_is_concrete():
    assert_concrete_property(ProblemType, "vars")

def test_elementwise_is_concrete():
    assert_concrete_property(ProblemType, "elementwise")

def test_elementwise_func_is_concrete():
    assert_concrete_property(ProblemType, "elementwise_func")

def test_elementwise_runner_is_concrete():
    assert_concrete_property(ProblemType, "elementwise_runner")

def test_replace_nan_values_by_is_concrete():
    assert_concrete_property(ProblemType, "replace_nan_values_by")

def test_exclude_from_serialization_is_concrete():
    assert_concrete_property(ProblemType, "exclude_from_serialization")

def test_callback_is_concrete():
    assert_concrete_property(ProblemType, "callback")

def test_strict_is_concrete():
    assert_concrete_property(ProblemType, "strict")

def test_data_is_concrete():
    assert_concrete_property(ProblemType, "data")

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(ProblemType, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ndecn_is_abstract():
    assert_abstract_property(ProblemType, "ndecn")

def test_decn_space_is_abstract():
    assert_abstract_property(ProblemType, "decn_space")

def test_decn_space_lower_is_abstract():
    assert_abstract_property(ProblemType, "decn_space_lower")

def test_decn_space_upper_is_abstract():
    assert_abstract_property(ProblemType, "decn_space_upper")

def test_nobj_is_abstract():
    assert_abstract_property(ProblemType, "nobj")

def test_obj_wt_is_abstract():
    assert_abstract_property(ProblemType, "obj_wt")

def test_nineqcv_is_abstract():
    assert_abstract_property(ProblemType, "nineqcv")

def test_ineqcv_wt_is_abstract():
    assert_abstract_property(ProblemType, "ineqcv_wt")

def test_neqcv_is_abstract():
    assert_abstract_property(ProblemType, "neqcv")

def test_eqcv_wt_is_abstract():
    assert_abstract_property(ProblemType, "eqcv_wt")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_evalfn_is_abstract(prob):
    assert_abstract_method(prob, "evalfn")

def test__evaluate_is_abstract(prob):
    assert_abstract_method(prob, "_evaluate")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_Problem_is_concrete():
    assert_concrete_function(check_is_ProblemType)

def test_check_is_Problem(prob):
    with not_raises(TypeError):
        check_is_ProblemType(prob, "prob")
    with pytest.raises(TypeError):
        check_is_ProblemType(None, "prob")
