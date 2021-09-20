import pytest
import numpy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import DenseTaxaMatrix
from pybropt.core.mat import is_DenseTaxaMatrix
from pybropt.core.mat import check_is_DenseTaxaMatrix
from pybropt.core.mat import cond_check_is_DenseTaxaMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([[3.3, 9.2, 5.6], [8.7, 3.7, 4.1], [9.0, 4.7, 3.8]])
    yield a

@pytest.fixture
def taxa_object():
    a = numpy.object_(["A", "B", "C"])
    yield a

@pytest.fixture
def taxa_grp_int64():
    a = numpy.int64([1,2,2])
    yield a

@pytest.fixture
def mat(mat_float64, taxa_object, taxa_grp_int64):
    yield DenseTaxaMatrix(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseTaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "__init__")

def test_adjoin_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "adjoin_taxa")

def test_delete_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "delete_taxa")

def test_insert_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "insert_taxa")

def test_select_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "select_taxa")

def test_concat_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "concat_taxa")

def test_append_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "append_taxa")

def test_remove_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "remove_taxa")

def test_incorp_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaMatrix, "incorp_taxa")

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_adjoin_taxa_cls(mat, mat_float64):
    m = mat.adjoin_taxa(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))

def test_adjoin_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_delete_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,2,None)
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_delete_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 0
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_delete_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,1,2]
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_insert_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,len(mat_float64),None)
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

def test_insert_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 1
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

def test_insert_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [e for e in range(len(mat_float64))]
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

def test_select_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))

def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [mat, mat]
    m = mat.concat_taxa(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
    assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

def test_append_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    mat.append_taxa(mat)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_append_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64):
    mat.append_taxa(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_remove_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,2,None)
    mat.remove_taxa(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_remove_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 0
    mat.remove_taxa(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_remove_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,1,2]
    mat.remove_taxa(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_incorp_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,len(mat_float64),None)
    mat.incorp_taxa(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

def test_incorp_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 1
    mat.incorp_taxa(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

def test_incorp_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [e for e in range(len(mat_float64))]
    mat.incorp_taxa(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_ntaxa_fget(mat, mat_float64):
    assert mat.ntaxa == len(mat_float64)

def test_ntaxa_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.ntaxa = len(mat_float64)

def test_ntaxa_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.ntaxa

def test_taxa_axis_fget(mat):
    assert mat.taxa_axis == 0

def test_taxa_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.taxa_axis = 1

def test_taxa_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_axis

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseTaxaMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseTaxaMatrix)

def test_is_DenseTaxaMatrix(mat):
    assert is_DenseTaxaMatrix(mat)

def test_check_is_DenseTaxaMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseTaxaMatrix)

def test_check_is_DenseTaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTaxaMatrix(None, "mat")

def test_cond_check_is_DenseTaxaMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseTaxaMatrix)
