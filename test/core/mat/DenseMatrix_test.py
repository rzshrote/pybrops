import os
from pathlib import Path
import pytest
import numpy
import h5py

from pybrops.test.assert_python import assert_class_documentation, assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DenseMatrix import DenseMatrix
from pybrops.core.mat.DenseMatrix import check_is_DenseMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.array([
        [3.3, 9.2, 5.6],
        [8.7, 3.7, 4.1],
        [9.0, 4.7, 3.8]
    ], dtype = float)
    yield a

@pytest.fixture
def mat(mat_float64):
    yield DenseMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseMatrix)

################################################################################
######################### Test special class functions #########################
################################################################################

### __init__

def test_init_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__init__")

####################################
###### Forward math operators ######

### __add__

def test___add___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__add__")

def test_add(mat, mat_float64):
    assert numpy.all((mat + 2.0) == (mat_float64 + 2.0))
    assert numpy.all((mat + mat_float64) == (mat_float64 + mat_float64))
    assert numpy.all((mat + mat) == (mat_float64 + mat_float64))

### __sub__

def test_sub_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__sub__")

def test_sub(mat, mat_float64):
    assert numpy.all((mat - 2.0) == (mat_float64 - 2.0))
    assert numpy.all((mat - mat_float64) == (mat_float64 - mat_float64))
    assert numpy.all((mat - mat) == (mat_float64 - mat_float64))

### __mul__

def test_mul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__mul__")

def test_mul(mat, mat_float64):
    assert numpy.all((mat * 2.0) == (mat_float64 * 2.0))
    assert numpy.all((mat * mat_float64) == (mat_float64 * mat_float64))
    assert numpy.all((mat * mat) == (mat_float64 * mat_float64))

### __matmul__

def test_matmul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__matmul__")

def test_matmul(mat, mat_float64):
    assert numpy.all((mat @ mat_float64) == (mat_float64 @ mat_float64))
    assert numpy.all((mat @ mat) == (mat_float64 @ mat_float64))

### __truediv__

def test_truediv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__truediv__")

def test_truediv(mat, mat_float64):
    assert numpy.all((mat / 2.0) == (mat_float64 / 2.0))
    assert numpy.all((mat / mat_float64) == (mat_float64 / mat_float64))
    assert numpy.all((mat / mat) == (mat_float64 / mat_float64))

### __floordiv__

def test_floordiv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__floordiv__")

def test_floordiv(mat, mat_float64):
    assert numpy.all((mat // 2.0) == (mat_float64 // 2.0))
    assert numpy.all((mat // mat_float64) == (mat_float64 // mat_float64))
    assert numpy.all((mat // mat) == (mat_float64 // mat_float64))

### __mod__

def test_mod_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__mod__")

### __divmod__

def test_divmod_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__divmod__")

### __pow__

def test_pow_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__pow__")

### __lshift__

def test_lshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__lshift__")

### __rshift__

def test_rshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rshift__")

### __and__

def test_and_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__and__")

### __xor__

def test_xor_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__xor__")

### __or__

def test_or_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__or__")

####################################
###### Reverse math operators ######

def test_radd_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__radd__")

def test_rsub_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rsub__")

def test_rmul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rmul__")

def test_rmatmul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rmatmul__")

def test_rtruediv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rtruediv__")

def test_rfloordiv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rfloordiv__")

def test_rmod_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rmod__")

def test_rdivmod_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rdivmod__")

def test_rlshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rlshift__")

def test_rrshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rrshift__")

def test_rand_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rand__")

def test_rxor_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__rxor__")

def test_ror_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ror__")

### Augment operators
def test_iadd_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__iadd__")

def test_isub_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__isub__")

def test_imul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__imul__")

def test_imatmul_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__imatmul__")

def test_itruediv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__itruediv__")

def test_ifloordiv_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ifloordiv__")

def test_imod_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__imod__")

def test_ipow_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ipow__")

def test_ilshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ilshift__")

def test_irshift_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__irshift__")

def test_iand_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__iand__")

def test_ixor_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ixor__")

def test_ior_is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ior__")

### Logical operators
def test___lt___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__lt__")

def test___le___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__le__")

def test___eq___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__eq__")

def test___ne___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ne__")

def test___gt___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__gt__")

def test___ge___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__ge__")

### Container operators
def test___len___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__len__")

def test___getitem___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__getitem__")

def test___setitem___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__setitem__")

def test___delitem___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__delitem__")

def test___iter___is_concrete():
    assert_method_isconcrete(DenseMatrix, "__iter__")

################################################################################
############################ Test Class Properties #############################
################################################################################
    
### mat

def test_mat_fget(mat, mat_float64):
    assert numpy.all(mat.mat == mat_float64)

################################################################################
############################# Test concrete methods ############################
################################################################################

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseMatrix, "to_hdf5")

def test_to_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_Path(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_h5py_File(mat):
    fp = "tmp.h5"
    h5file = h5py.File(fp, "a")
    with not_raises(Exception):
        mat.to_hdf5(h5file)
    h5file.close()
    assert os.path.exists(fp)
    os.remove(fp)

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

### from_hdf5

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseMatrix, "from_hdf5")

def test_from_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    out = DenseMatrix.from_hdf5(fp)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    os.remove(fp)

def test_from_hdf5_Path(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    out = DenseMatrix.from_hdf5(fp)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    os.remove(fp)

def test_from_hdf5_h5py_File(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseMatrix.from_hdf5(h5file)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    h5file.close()
    os.remove(fp)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseMatrix)

def test_check_is_DenseMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMatrix(None, "mat")
