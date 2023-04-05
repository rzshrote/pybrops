import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.prob.Problem import Problem, check_is_Problem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield Problem()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(Problem)

################################################################################
########################### Test concrete properties ###########################
################################################################################
def test_n_var_is_concrete():
    assert_concrete_property(Problem, "n_var")

def test_n_obj_is_concrete():
    assert_concrete_property(Problem, "n_obj")

def test_n_ieq_constr_is_concrete():
    assert_concrete_property(Problem, "n_ieq_constr")

def test_n_eq_constr_is_concrete():
    assert_concrete_property(Problem, "n_eq_constr")

def test_xl_is_concrete():
    assert_concrete_property(Problem, "xl")

def test_xu_is_concrete():
    assert_concrete_property(Problem, "xu")

def test_vtype_is_concrete():
    assert_concrete_property(Problem, "vtype")

def test_vars_is_concrete():
    assert_concrete_property(Problem, "vars")

def test_elementwise_is_concrete():
    assert_concrete_property(Problem, "elementwise")

def test_elementwise_func_is_concrete():
    assert_concrete_property(Problem, "elementwise_func")

def test_elementwise_runner_is_concrete():
    assert_concrete_property(Problem, "elementwise_runner")

def test_replace_nan_values_by_is_concrete():
    assert_concrete_property(Problem, "replace_nan_values_by")

def test_exclude_from_serialization_is_concrete():
    assert_concrete_property(Problem, "exclude_from_serialization")

def test_callback_is_concrete():
    assert_concrete_property(Problem, "callback")

def test_strict_is_concrete():
    assert_concrete_property(Problem, "strict")

def test_data_is_concrete():
    assert_concrete_property(Problem, "data")

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(Problem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ndecn_is_abstract():
    assert_abstract_property(Problem, "ndecn")

def test_decn_space_is_abstract():
    assert_abstract_property(Problem, "decn_space")

def test_decn_space_lower_is_abstract():
    assert_abstract_property(Problem, "decn_space_lower")

def test_decn_space_upper_is_abstract():
    assert_abstract_property(Problem, "decn_space_upper")

def test_nobj_is_abstract():
    assert_abstract_property(Problem, "nobj")

def test_obj_wt_is_abstract():
    assert_abstract_property(Problem, "obj_wt")

def test_nineqcv_is_abstract():
    assert_abstract_property(Problem, "nineqcv")

def test_ineqcv_wt_is_abstract():
    assert_abstract_property(Problem, "ineqcv_wt")

def test_neqcv_is_abstract():
    assert_abstract_property(Problem, "neqcv")

def test_eqcv_wt_is_abstract():
    assert_abstract_property(Problem, "eqcv_wt")

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
    assert_concrete_function(check_is_Problem)

def test_check_is_Problem(prob):
    with not_raises(TypeError):
        check_is_Problem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_Problem(None, "prob")
