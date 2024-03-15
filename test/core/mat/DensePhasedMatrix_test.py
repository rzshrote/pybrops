import pytest
import numpy

from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DensePhasedMatrix import DensePhasedMatrix
from pybrops.core.mat.DensePhasedMatrix import check_is_DensePhasedMatrix

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
    assert_class_documentation(DensePhasedMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "__init__")

def test_adjoin_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "adjoin_phase")

def test_delete_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "delete_phase")

def test_insert_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "insert_phase")

def test_select_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "select_phase")

def test_concat_phase_is_concrete():
    assert_classmethod_isconcrete(DensePhasedMatrix, "concat_phase")

def test_append_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "append_phase")

def test_remove_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "remove_phase")

def test_incorp_phase_is_concrete():
    assert_method_isconcrete(DensePhasedMatrix, "incorp_phase")

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
def test_check_is_DensePhasedMatrix_is_concrete():
    assert_function_isconcrete(check_is_DensePhasedMatrix)

def test_check_is_DensePhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_DensePhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DensePhasedMatrix(None, "mat")
