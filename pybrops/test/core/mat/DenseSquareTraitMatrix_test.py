import pytest
import numpy
import copy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.DenseSquareTraitMatrix import DenseSquareTraitMatrix
from pybrops.core.mat.DenseSquareTraitMatrix import check_is_DenseSquareTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntrait():
    yield 4

@pytest.fixture
def mat_float64(ntrait):
    out = numpy.random.random((ntrait,ntrait))
    yield out

@pytest.fixture
def trait_object(ntrait):
    out = numpy.array(["Trait"+str(i) for i in range(ntrait)], dtype = object)
    yield out

@pytest.fixture
def trait_lexsort_indices(trait_object):
    out = numpy.lexsort((trait_object,))
    yield out

@pytest.fixture
def mat(mat_float64, trait_object):
    out = DenseSquareTraitMatrix(
        mat = mat_float64,
        trait = trait_object,
    )
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseSquareTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "__init__")

def test_adjoin_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "adjoin")

def test_adjoin_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "adjoin_trait")

def test_delete_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "delete")

def test_delete_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "delete_trait")

def test_insert_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "insert")

def test_insert_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "insert_trait")

def test_select_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "select")

def test_select_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "select_trait")

def test_concat_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "concat")

def test_concat_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "concat_trait")

def test_append_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "append")

def test_append_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "append_trait")

def test_remove_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "remove")

def test_remove_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "remove_trait")

def test_incorp_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "incorp")

def test_incorp_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "incorp_trait")

def test_lexsort_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "lexsort")

def test_lexsort_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "lexsort_trait")

def test_reorder_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "reorder")

def test_reorder_trait_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "reorder_trait")

def test_sort_is_concrete():
    assert_concrete_method(DenseSquareTraitMatrix, "sort")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########### Matrix element copy-on-manipulation ############
def test_adjoin_trait_cls(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    mtrue = numpy.block([[a,b],[b,a]])
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_adjoin_trait_ndarray(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat_float64, trait = trait_object)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    mtrue = numpy.block([[a,b],[b,a]])
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_delete_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,2,None)
    m = mat.delete_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_trait_cls_int(mat, mat_float64, trait_object):
    obj = 0
    m = mat.delete_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,1,2]
    m = mat.delete_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

# def test_insert_trait_cls_slice(mat, mat_float64, trait_object):
#     obj = slice(0,len(mat_float64),None)
#     m = mat.insert_trait(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
#
# def test_insert_trait_cls_int(mat, mat_float64, trait_object):
#     obj = 1
#     m = mat.insert_trait(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
#
# def test_insert_trait_cls_array_like(mat, mat_float64, trait_object):
#     obj = [e for e in range(len(mat_float64))]
#     m = mat.insert_trait(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

def test_select_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,0,1]
    m = mat.select_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.take(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))

# def test_concat_trait_cls(mat, mat_float64, trait_object):
#     obj = [mat, mat]
#     m = mat.concat_trait(obj)
#     assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.trait_axis))
#     assert numpy.all(m.trait == numpy.concatenate([trait_object, trait_object], axis = 0))

########### Matrix element in-place-manipulation ###########
def test_append_trait_cls(mat, mat_float64, trait_object):
    mat.append_trait(mat)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    mtrue = numpy.block([[a,b],[b,a]])
    assert numpy.array_equal(mat.mat, mtrue, equal_nan = True)
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_append_trait_ndarray(mat, mat_float64, trait_object):
    mat.append_trait(mat_float64, trait = trait_object)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    mtrue = numpy.block([[a,b],[b,a]])
    assert numpy.array_equal(mat.mat, mtrue, equal_nan = True)
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_remove_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,2,None)
    mat.remove_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_trait_cls_int(mat, mat_float64, trait_object):
    obj = 0
    mat.remove_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,1,2]
    mat.remove_trait(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.trait == numpy.delete(trait_object, obj, axis = 0))

# def test_incorp_trait_cls_slice(mat, mat_float64, trait_object):
#     obj = slice(0,len(mat_float64),None)
#     mat.incorp_trait(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
#
# def test_incorp_trait_cls_int(mat, mat_float64, trait_object):
#     obj = 1
#     mat.incorp_trait(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
#
# def test_incorp_trait_cls_array_like(mat, mat_float64, trait_object):
#     obj = [e for e in range(len(mat_float64))]
#     mat.incorp_trait(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
#     assert numpy.all(mat.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))

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
def test_check_is_DenseSquareTraitMatrix_is_concrete():
    assert_concrete_function(check_is_DenseSquareTraitMatrix)

def test_check_is_DenseSquareTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseSquareTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseSquareTraitMatrix(None, "mat")
