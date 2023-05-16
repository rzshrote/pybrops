import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.DenseMatrix import DenseMatrix
from pybrops.core.mat.DenseMatrix import is_DenseMatrix
from pybrops.core.mat.DenseMatrix import check_is_DenseMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6],
        [8.7, 3.7, 4.1],
        [9.0, 4.7, 3.8]
    ])
    yield a

@pytest.fixture
def mat(mat_float64):
    yield DenseMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseMatrix, "__init__")

### Forward math operators
def test_add_is_concrete():
    assert_concrete_method(DenseMatrix, "__add__")

def test_sub_is_concrete():
    assert_concrete_method(DenseMatrix, "__sub__")

def test_mul_is_concrete():
    assert_concrete_method(DenseMatrix, "__mul__")

def test_matmul_is_concrete():
    assert_concrete_method(DenseMatrix, "__matmul__")

def test_truediv_is_concrete():
    assert_concrete_method(DenseMatrix, "__truediv__")

def test_floordiv_is_concrete():
    assert_concrete_method(DenseMatrix, "__floordiv__")

def test_mod_is_concrete():
    assert_concrete_method(DenseMatrix, "__mod__")

def test_divmod_is_concrete():
    assert_concrete_method(DenseMatrix, "__divmod__")

def test_pow_is_concrete():
    assert_concrete_method(DenseMatrix, "__pow__")

def test_lshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__lshift__")

def test_rshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__rshift__")

def test_and_is_concrete():
    assert_concrete_method(DenseMatrix, "__and__")

def test_xor_is_concrete():
    assert_concrete_method(DenseMatrix, "__xor__")

def test_or_is_concrete():
    assert_concrete_method(DenseMatrix, "__or__")

### Reverse math operators
def test_radd_is_concrete():
    assert_concrete_method(DenseMatrix, "__radd__")

def test_rsub_is_concrete():
    assert_concrete_method(DenseMatrix, "__rsub__")

def test_rmul_is_concrete():
    assert_concrete_method(DenseMatrix, "__rmul__")

def test_rmatmul_is_concrete():
    assert_concrete_method(DenseMatrix, "__rmatmul__")

def test_rtruediv_is_concrete():
    assert_concrete_method(DenseMatrix, "__rtruediv__")

def test_rfloordiv_is_concrete():
    assert_concrete_method(DenseMatrix, "__rfloordiv__")

def test_rmod_is_concrete():
    assert_concrete_method(DenseMatrix, "__rmod__")

def test_rdivmod_is_concrete():
    assert_concrete_method(DenseMatrix, "__rdivmod__")

def test_rlshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__rlshift__")

def test_rrshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__rrshift__")

def test_rand_is_concrete():
    assert_concrete_method(DenseMatrix, "__rand__")

def test_rxor_is_concrete():
    assert_concrete_method(DenseMatrix, "__rxor__")

def test_ror_is_concrete():
    assert_concrete_method(DenseMatrix, "__ror__")

### Augment operators
def test_iadd_is_concrete():
    assert_concrete_method(DenseMatrix, "__iadd__")

def test_isub_is_concrete():
    assert_concrete_method(DenseMatrix, "__isub__")

def test_imul_is_concrete():
    assert_concrete_method(DenseMatrix, "__imul__")

def test_imatmul_is_concrete():
    assert_concrete_method(DenseMatrix, "__imatmul__")

def test_itruediv_is_concrete():
    assert_concrete_method(DenseMatrix, "__itruediv__")

def test_ifloordiv_is_concrete():
    assert_concrete_method(DenseMatrix, "__ifloordiv__")

def test_imod_is_concrete():
    assert_concrete_method(DenseMatrix, "__imod__")

def test_ipow_is_concrete():
    assert_concrete_method(DenseMatrix, "__ipow__")

def test_ilshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__ilshift__")

def test_irshift_is_concrete():
    assert_concrete_method(DenseMatrix, "__irshift__")

def test_iand_is_concrete():
    assert_concrete_method(DenseMatrix, "__iand__")

def test_ixor_is_concrete():
    assert_concrete_method(DenseMatrix, "__ixor__")

def test_ior_is_concrete():
    assert_concrete_method(DenseMatrix, "__ior__")


################################################################################
############################ Test Class Properties #############################
################################################################################
def test_mat_fget(mat, mat_float64):
    assert numpy.all(mat.mat == mat_float64)

################################################################################
######################### Test special class functions #########################
################################################################################
def test_add(mat, mat_float64):
    assert numpy.all((mat + 2.0) == (mat_float64 + 2.0))
    assert numpy.all((mat + mat_float64) == (mat_float64 + mat_float64))
    assert numpy.all((mat + mat) == (mat_float64 + mat_float64))

def test_sub(mat, mat_float64):
    assert numpy.all((mat - 2.0) == (mat_float64 - 2.0))
    assert numpy.all((mat - mat_float64) == (mat_float64 - mat_float64))
    assert numpy.all((mat - mat) == (mat_float64 - mat_float64))

def test_mul(mat, mat_float64):
    assert numpy.all((mat * 2.0) == (mat_float64 * 2.0))
    assert numpy.all((mat * mat_float64) == (mat_float64 * mat_float64))
    assert numpy.all((mat * mat) == (mat_float64 * mat_float64))

def test_matmul(mat, mat_float64):
    assert numpy.all((mat @ mat_float64) == (mat_float64 @ mat_float64))
    assert numpy.all((mat @ mat) == (mat_float64 @ mat_float64))

def test_truediv(mat, mat_float64):
    assert numpy.all((mat / 2.0) == (mat_float64 / 2.0))
    assert numpy.all((mat / mat_float64) == (mat_float64 / mat_float64))
    assert numpy.all((mat / mat) == (mat_float64 / mat_float64))

def test_floordiv(mat, mat_float64):
    assert numpy.all((mat // 2.0) == (mat_float64 // 2.0))
    assert numpy.all((mat // mat_float64) == (mat_float64 // mat_float64))
    assert numpy.all((mat // mat) == (mat_float64 // mat_float64))

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseMatrix_is_concrete():
    assert_concrete_function(is_DenseMatrix)

def test_is_DenseMatrix(mat):
    assert is_DenseMatrix(mat)

def test_check_is_DenseMatrix_is_concrete():
    assert_concrete_function(check_is_DenseMatrix)

def test_check_is_DenseMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMatrix(None, "mat")
