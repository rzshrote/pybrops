import pytest
import numpy
import copy

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.core.mat.DenseTraitMatrix import is_DenseTraitMatrix
from pybrops.core.mat.DenseTraitMatrix import check_is_DenseTraitMatrix
from pybrops.core.mat.DenseTraitMatrix import cond_check_is_DenseTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([[3.3, 9.2, 5.6], [8.7, 3.7, 4.1], [9.0, 4.7, 3.8]])
    yield a

@pytest.fixture
def trait_object():
    a = numpy.object_(["yield", "oil", "protein"])
    yield a

@pytest.fixture
def trait_lexsort_indices(trait_object):
    a = numpy.lexsort((trait_object,))
    yield a

@pytest.fixture
def mat(mat_float64, trait_object):
    out = DenseTraitMatrix(mat_float64, trait = trait_object)
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "__init__")

def test_copy_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "__copy__")

def test_deepcopy_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "__deepcopy__")

def test_adjoin_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "adjoin_trait")

def test_delete_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "delete_trait")

def test_insert_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "insert_trait")

def test_select_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "select_trait")

def test_concat_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "concat_trait")

def test_append_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "append_trait")

def test_remove_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "remove_trait")

def test_incorp_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "incorp_trait")

def test_lexsort_trait_is_concrete():
    generic_assert_concrete_method(DenseTraitMatrix, "lexsort_trait")

# TODO: # FIXME: not_raises fails for an edge case
# def test_sort_trait_is_concrete():
#     generic_assert_abstract_method(DenseTraitMatrix, "sort_trait")
#
# def test_group_trait_is_concrete():
#     generic_assert_abstract_method(DenseTraitMatrix, "group_trait")
#
# def test_is_grouped_trait_is_concrete():
#     generic_assert_abstract_method(DenseTraitMatrix, "is_grouped_trait")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test_copy(mat):
    m = copy.copy(mat)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.trait == mat.trait)

def test_deepcopy(mat):
    m = copy.deepcopy(mat)
    # make sure object ID's are different
    assert id(m.mat) != id(mat.mat)
    assert id(m.trait) != id(mat.trait)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.trait == mat.trait)

################################################################################
############################ Test Class Properties #############################
################################################################################

################### Taxa Data Properites ###################
def test_trait_fget(mat, trait_object):
    assert numpy.all(mat.trait == trait_object)

def test_trait_fset(mat, trait_object):
    mat.trait = trait_object
    assert numpy.all(mat.trait == trait_object)

def test_trait_fdel(mat):
    del mat.trait
    with pytest.raises(AttributeError):
        mat.trait

################# Taxa Metadata Properites #################
def test_ntrait_fget(mat, mat_float64):
    assert mat.ntrait == len(mat_float64)

def test_ntrait_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.ntrait = len(mat_float64)

def test_ntrait_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.ntrait

def test_trait_axis_fget(mat):
    assert mat.trait_axis == 0

def test_trait_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.trait_axis = 1

def test_trait_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.trait_axis

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########### Matrix element copy-on-manipulation ############
def test_adjoin_trait_cls(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_adjoin_trait_ndarray(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat_float64, trait = trait_object)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_delete_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,2,None)
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_trait_cls_int(mat, mat_float64, trait_object):
    obj = 0
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,1,2]
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_insert_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,len(mat_float64),None)
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_insert_trait_cls_int(mat, mat_float64, trait_object):
    obj = 1
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_insert_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [e for e in range(len(mat_float64))]
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_select_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,0,1]
    m = mat.select_trait(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))

def test_concat_trait_cls(mat, mat_float64, trait_object):
    obj = [mat, mat]
    m = mat.concat_trait(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.concatenate([trait_object, trait_object], axis = 0))

########### Matrix element in-place-manipulation ###########
def test_append_trait_cls(mat, mat_float64, trait_object):
    mat.append_trait(mat)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_append_trait_ndarray(mat, mat_float64, trait_object):
    mat.append_trait(mat_float64, trait = trait_object)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_remove_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,2,None)
    mat.remove_trait(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_trait_cls_int(mat, mat_float64, trait_object):
    obj = 0
    mat.remove_trait(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,1,2]
    mat.remove_trait(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

def test_incorp_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,len(mat_float64),None)
    mat.incorp_trait(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_incorp_trait_cls_int(mat, mat_float64, trait_object):
    obj = 1
    mat.incorp_trait(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_incorp_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [e for e in range(len(mat_float64))]
    mat.incorp_trait(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

##################### Sorting Methods ######################
def test_lexsort_trait_None(mat, trait_lexsort_indices):
    ix = mat.lexsort_trait(keys = None)
    assert numpy.all(ix == trait_lexsort_indices)

def test_lexsort_trait_tuple(mat, trait_object, trait_lexsort_indices):
    ix = mat.lexsort_trait(keys = (trait_object,))
    assert numpy.all(ix == trait_lexsort_indices)

def test_reorder_trait_array_like(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.reorder_trait(trait_lexsort_indices)
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_sort_trait_None(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.sort_trait(keys = None)
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_sort_trait_tuple(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.sort_trait(keys = (trait_object,))
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseTraitMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseTraitMatrix)

def test_is_DenseTraitMatrix(mat):
    assert is_DenseTraitMatrix(mat)

def test_check_is_DenseTraitMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseTraitMatrix)

def test_check_is_DenseTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTraitMatrix(None, "mat")

def test_cond_check_is_DenseTraitMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseTraitMatrix)
