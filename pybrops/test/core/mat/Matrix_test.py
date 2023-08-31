import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.Matrix import check_is_Matrix

from pybrops.test.core.mat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(Matrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(Matrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_mat_is_abstract():
    assert_abstract_property(Matrix, "mat")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_Matrix_add_is_abstract():
    assert_abstract_method(Matrix, "__add__")

def test_Matrix_sub_is_abstract():
    assert_abstract_method(Matrix, "__sub__")

def test_Matrix_mul_is_abstract():
    assert_abstract_method(Matrix, "__mul__")

def test_Matrix_matmul_is_abstract():
    assert_abstract_method(Matrix, "__matmul__")

def test_Matrix_truediv_is_abstract():
    assert_abstract_method(Matrix, "__truediv__")

def test_Matrix_floordiv_is_abstract():
    assert_abstract_method(Matrix, "__floordiv__")

def test_Matrix_mod_is_abstract():
    assert_abstract_method(Matrix, "__mod__")

def test_Matrix_divmod_is_abstract():
    assert_abstract_method(Matrix, "__divmod__")

def test_Matrix_pow_is_abstract():
    assert_abstract_method(Matrix, "__pow__")

def test_Matrix_lshift_is_abstract():
    assert_abstract_method(Matrix, "__lshift__")

def test_Matrix_rshift_is_abstract():
    assert_abstract_method(Matrix, "__rshift__")

def test_Matrix_and_is_abstract():
    assert_abstract_method(Matrix, "__and__")

def test_Matrix_xor_is_abstract():
    assert_abstract_method(Matrix, "__xor__")

def test_Matrix_or_is_abstract():
    assert_abstract_method(Matrix, "__or__")

def test_Matrix_rand_is_abstract():
    assert_abstract_method(Matrix, "__rand__")

def test_Matrix_rsub_is_abstract():
    assert_abstract_method(Matrix, "__rsub__")

def test_Matrix_rmul_is_abstract():
    assert_abstract_method(Matrix, "__rmul__")

def test_Matrix_rmatmul_is_abstract():
    assert_abstract_method(Matrix, "__rmatmul__")

def test_Matrix_rtruediv_is_abstract():
    assert_abstract_method(Matrix, "__rtruediv__")

def test_Matrix_rfloordiv_is_abstract():
    assert_abstract_method(Matrix, "__rfloordiv__")

def test_Matrix_rmod_is_abstract():
    assert_abstract_method(Matrix, "__rmod__")

def test_Matrix_rlshift_is_abstract():
    assert_abstract_method(Matrix, "__rlshift__")

def test_Matrix_rrshift_is_abstract():
    assert_abstract_method(Matrix, "__rrshift__")

def test_Matrix_rand_is_abstract():
    assert_abstract_method(Matrix, "__rand__")

def test_Matrix_rxor_is_abstract():
    assert_abstract_method(Matrix, "__rxor__")

def test_Matrix_ror_is_abstract():
    assert_abstract_method(Matrix, "__ror__")

def test_Matrix_iadd_is_abstract():
    assert_abstract_method(Matrix, "__iadd__")

def test_Matrix_isub_is_abstract():
    assert_abstract_method(Matrix, "__isub__")

def test_Matrix_imul_is_abstract():
    assert_abstract_method(Matrix, "__imul__")

def test_Matrix_imatmul_is_abstract():
    assert_abstract_method(Matrix, "__imatmul__")

def test_Matrix_itruediv_is_abstract():
    assert_abstract_method(Matrix, "__itruediv__")

def test_Matrix_ifloordiv_is_abstract():
    assert_abstract_method(Matrix, "__ifloordiv__")

def test_Matrix_imod_is_abstract():
    assert_abstract_method(Matrix, "__imod__")

def test_Matrix_ipow_is_abstract():
    assert_abstract_method(Matrix, "__ipow__")

def test_Matrix_ilshift_is_abstract():
    assert_abstract_method(Matrix, "__ilshift__")

def test_Matrix_irshift_is_abstract():
    assert_abstract_method(Matrix, "__irshift__")

def test_Matrix_iand_is_abstract():
    assert_abstract_method(Matrix, "__iand__")

def test_Matrix_ixor_is_abstract():
    assert_abstract_method(Matrix, "__ixor__")

def test_Matrix_ior_is_abstract():
    assert_abstract_method(Matrix, "__ior__")

def test_Matrix_lt_is_abstract():
    assert_abstract_method(Matrix, "__lt__")

def test_Matrix_le_is_abstract():
    assert_abstract_method(Matrix, "__le__")

def test_Matrix_eq_is_abstract():
    assert_abstract_method(Matrix, "__eq__")

def test_Matrix_ne_is_abstract():
    assert_abstract_method(Matrix, "__ne__")

def test_Matrix_gt_is_abstract():
    assert_abstract_method(Matrix, "__gt__")

def test_Matrix_ge_is_abstract():
    assert_abstract_method(Matrix, "__ge__")

def test_Matrix_len_is_abstract():
    assert_abstract_method(Matrix, "__len__")

def test_Matrix_getitem_is_abstract():
    assert_abstract_method(Matrix, "__getitem__")

def test_Matrix_delitem_is_abstract():
    assert_abstract_method(Matrix, "__delitem__")

def test_Matrix_iter_is_abstract():
    assert_abstract_method(Matrix, "__iter__")

def test_Matrix_copy_is_abstract():
    assert_abstract_method(Matrix, "__copy__")
    assert_abstract_method(Matrix, "copy")

def test_Matrix_deepcopy_is_abstract():
    assert_abstract_method(Matrix, "__deepcopy__")
    assert_abstract_method(Matrix, "deepcopy")

def test_Matrix_adjoin_is_abstract():
    assert_abstract_method(Matrix, "adjoin")

def test_Matrix_delete_is_abstract():
    assert_abstract_method(Matrix, "delete")

def test_Matrix_insert_is_abstract():
    assert_abstract_method(Matrix, "insert")

def test_Matrix_select_is_abstract():
    assert_abstract_method(Matrix, "select")

def test_Matrix_concat_is_abstract():
    assert_abstract_method(Matrix, "concat")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_Matrix_is_concrete():
    assert_concrete_function(check_is_Matrix)

def test_check_is_Matrix(mat):
    with not_raises(TypeError):
        check_is_Matrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_Matrix(None, "mat")
