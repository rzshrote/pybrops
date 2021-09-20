import pytest
import numpy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import DensePhasedMatrix
from pybropt.core.mat import is_DensePhasedMatrix
from pybropt.core.mat import check_is_DensePhasedMatrix
from pybropt.core.mat import cond_check_is_DensePhasedMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([[3.3, 9.2, 5.6], [8.7, 3.7, 4.1], [9.0, 4.7, 3.8]])
    yield a

@pytest.fixture
def mat(mat_float64):
    yield DensePhasedMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DensePhasedMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "__init__")

def test_adjoin_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "adjoin_phase")

def test_delete_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "delete_phase")

def test_insert_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "insert_phase")

def test_select_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "select_phase")

def test_concat_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "concat_phase")

def test_append_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "append_phase")

def test_remove_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "remove_phase")

def test_incorp_phase_is_concrete():
    generic_assert_concrete_method(DensePhasedMatrix, "incorp_phase")

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_adjoin_phase_cls(mat, mat_float64):
    m = mat.adjoin_phase(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.phase_axis))

def test_adjoin_phase_ndarray(mat, mat_float64):
    m = mat.adjoin_phase(mat_float64)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.phase_axis))

def test_delete_phase_cls_slice(mat, mat_float64):
    obj = slice(0,2,None)
    m = mat.delete_phase(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_delete_phase_cls_int(mat, mat_float64):
    obj = 0
    m = mat.delete_phase(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_delete_phase_cls_array_like(mat, mat_float64):
    obj = [0,1,2]
    m = mat.delete_phase(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_insert_phase_cls_slice(mat, mat_float64):
    obj = slice(0,len(mat_float64),None)
    m = mat.insert_phase(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

def test_insert_phase_cls_int(mat, mat_float64):
    obj = 1
    m = mat.insert_phase(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

def test_insert_phase_cls_array_like(mat, mat_float64):
    obj = [e for e in range(len(mat_float64))]
    m = mat.insert_phase(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

def test_select_phase_cls_array_like(mat, mat_float64):
    obj = [0,0,1]
    m = mat.select_phase(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.phase_axis))

def test_concat_phase_cls(mat, mat_float64):
    obj = [mat, mat]
    m = mat.concat_phase(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.phase_axis))

def test_append_phase_cls(mat, mat_float64):
    mat.append_phase(mat)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.phase_axis))

def test_append_phase_ndarray(mat, mat_float64):
    mat.append_phase(mat_float64)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.phase_axis))

def test_remove_phase_cls_slice(mat, mat_float64):
    obj = slice(0,2,None)
    mat.remove_phase(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_remove_phase_cls_int(mat, mat_float64):
    obj = 0
    mat.remove_phase(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_remove_phase_cls_array_like(mat, mat_float64):
    obj = [0,1,2]
    mat.remove_phase(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.phase_axis))

def test_incorp_phase_cls_slice(mat, mat_float64):
    obj = slice(0,len(mat_float64),None)
    mat.incorp_phase(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

def test_incorp_phase_cls_int(mat, mat_float64):
    obj = 1
    mat.incorp_phase(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

def test_incorp_phase_cls_array_like(mat, mat_float64):
    obj = [e for e in range(len(mat_float64))]
    mat.incorp_phase(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.phase_axis))

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_nphase_fget(mat, mat_float64):
    assert mat.nphase == len(mat_float64)

def test_nphase_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.nphase = len(mat_float64)

def test_nphase_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.nphase

def test_phase_axis_fget(mat):
    assert mat.phase_axis == 0

def test_phase_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.phase_axis = 1

def test_phase_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.phase_axis

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DensePhasedMatrix_is_concrete():
    generic_assert_concrete_function(is_DensePhasedMatrix)

def test_is_DensePhasedMatrix(mat):
    assert is_DensePhasedMatrix(mat)

def test_check_is_DensePhasedMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DensePhasedMatrix)

def test_check_is_DensePhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_DensePhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DensePhasedMatrix(None, "mat")

def test_cond_check_is_DensePhasedMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DensePhasedMatrix)
