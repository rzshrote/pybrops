import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import Matrix
from pybropt.core.mat import is_Matrix
from pybropt.core.mat import check_is_Matrix
from pybropt.core.mat import cond_check_is_Matrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield Matrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(Matrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(Matrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_add_is_abstract(mat):
    generic_assert_abstract_method(mat, "__add__")

def test_sub_is_abstract(mat):
    generic_assert_abstract_method(mat, "__sub__")

def test_mul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__mul__")

def test_matmul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__matmul__")

def test_truediv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__truediv__")

def test_floordiv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__floordiv__")

def test_mod_is_abstract(mat):
    generic_assert_abstract_method(mat, "__mod__")

def test_divmod_is_abstract(mat):
    generic_assert_abstract_method(mat, "__divmod__")

def test_pow_is_abstract(mat):
    generic_assert_abstract_method(mat, "__pow__")

def test_lshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__lshift__")

def test_rshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rshift__")

def test_and_is_abstract(mat):
    generic_assert_abstract_method(mat, "__and__")

def test_xor_is_abstract(mat):
    generic_assert_abstract_method(mat, "__xor__")

def test_or_is_abstract(mat):
    generic_assert_abstract_method(mat, "__or__")

def test_rand_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rand__")

def test_rsub_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rsub__")

def test_rmul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rmul__")

def test_rmatmul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rmatmul__")

def test_rtruediv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rtruediv__")

def test_rfloordiv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rfloordiv__")

def test_rmod_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rmod__")

def test_rlshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rlshift__")

def test_rrshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rrshift__")

def test_rand_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rand__")

def test_rxor_is_abstract(mat):
    generic_assert_abstract_method(mat, "__rxor__")

def test_ror_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ror__")

def test_iadd_is_abstract(mat):
    generic_assert_abstract_method(mat, "__iadd__")

def test_isub_is_abstract(mat):
    generic_assert_abstract_method(mat, "__isub__")

def test_imul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__imul__")

def test_imatmul_is_abstract(mat):
    generic_assert_abstract_method(mat, "__imatmul__")

def test_itruediv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__itruediv__")

def test_ifloordiv_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ifloordiv__")

def test_imod_is_abstract(mat):
    generic_assert_abstract_method(mat, "__imod__")

def test_ipow_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ipow__")

def test_ilshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ilshift__")

def test_irshift_is_abstract(mat):
    generic_assert_abstract_method(mat, "__irshift__")

def test_iand_is_abstract(mat):
    generic_assert_abstract_method(mat, "__iand__")

def test_ixor_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ixor__")

def test_ior_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ior__")

def test_lt_is_abstract(mat):
    generic_assert_abstract_method(mat, "__lt__")

def test_le_is_abstract(mat):
    generic_assert_abstract_method(mat, "__le__")

def test_eq_is_abstract(mat):
    generic_assert_abstract_method(mat, "__eq__")

def test_ne_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ne__")

def test_gt_is_abstract(mat):
    generic_assert_abstract_method(mat, "__gt__")

def test_ge_is_abstract(mat):
    generic_assert_abstract_method(mat, "__ge__")

def test_len_is_abstract(mat):
    generic_assert_abstract_method(mat, "__len__")

def test_getitem_is_abstract(mat):
    generic_assert_abstract_method(mat, "__getitem__")

def test_delitem_is_abstract(mat):
    generic_assert_abstract_method(mat, "__delitem__")

def test_iter_is_abstract(mat):
    generic_assert_abstract_method(mat, "__iter__")

def test_copy_is_abstract(mat):
    generic_assert_abstract_method(mat, "__copy__")
    generic_assert_abstract_method(mat, "copy")

def test_deepcopy_is_abstract(mat):
    generic_assert_abstract_method(mat, "__deepcopy__")
    generic_assert_abstract_method(mat, "deepcopy")

def test_mat_is_abstract():
    generic_assert_abstract_property(Matrix, "mat")

def test_adjoin_is_abstract(mat):
    generic_assert_abstract_method(mat, "adjoin")

def test_delete_is_abstract(mat):
    generic_assert_abstract_method(mat, "delete")

def test_insert_is_abstract(mat):
    generic_assert_abstract_method(mat, "insert")

def test_select_is_abstract(mat):
    generic_assert_abstract_method(mat, "select")

def test_concat_is_abstract(mat):
    generic_assert_abstract_method(mat, "concat")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_Matrix_is_concrete():
    generic_assert_concrete_function(is_Matrix)

def test_is_Matrix(mat):
    assert is_Matrix(mat)

def test_check_is_Matrix_is_concrete():
    generic_assert_concrete_function(check_is_Matrix)

def test_check_is_Matrix(mat):
    with not_raises(TypeError):
        check_is_Matrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_Matrix(None, "mat")

def test_cond_check_is_Matrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_Matrix)

# def test_cond_check_is_Matrix():
#     with not_raises(TypeError):
#         cond_check_is_Matrix(mat, "mat")
#     # with not_raises(TypeError):
#     #     cond_check_is_Matrix(None, "mat")
#     with pytest.raises(TypeError):
#         cond_check_is_Matrix("string", "mat")
