import os
from pathlib import Path
import pytest
import numpy
import copy
import h5py

from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.core.mat.DenseTraitMatrix import check_is_DenseTraitMatrix

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
    assert_class_documentation(DenseTraitMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "__init__")

### __copy__

def test___copy___is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "__deepcopy__")

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
    with pytest.raises(AttributeError):
        del mat.trait

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
############################# Test concrete methods ############################
################################################################################

### copy

def test_copy_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "copy")

def test_copy(mat):
    m = copy.copy(mat)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.trait == mat.trait)

### deepcopy

def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "deepcopy")

def test_deepcopy(mat):
    m = copy.deepcopy(mat)
    # make sure object ID's are different
    assert id(m.mat) != id(mat.mat)
    assert id(m.trait) != id(mat.trait)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.trait == mat.trait)

############################################################
########### Matrix element copy-on-manipulation ############

### adjoin_trait

def test_adjoin_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "adjoin_trait")

def test_adjoin_trait_cls(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_adjoin_trait_ndarray(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat_float64, trait = trait_object)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

### delete_trait

def test_delete_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "delete_trait")

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

### insert_trait

def test_insert_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "insert_trait")

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

### select_trait

def test_select_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "select_trait")

def test_select_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,0,1]
    m = mat.select_trait(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))

### concat_trait

def test_concat_trait_is_concrete():
    assert_classmethod_isconcrete(DenseTraitMatrix, "concat_trait")

def test_concat_trait_cls(mat, mat_float64, trait_object):
    obj = [mat, mat]
    m = mat.concat_trait(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.concatenate([trait_object, trait_object], axis = 0))

############################################################
########### Matrix element in-place-manipulation ###########

### append_trait

def test_append_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "append_trait")

def test_append_trait_cls(mat, mat_float64, trait_object):
    mat.append_trait(mat)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_append_trait_ndarray(mat, mat_float64, trait_object):
    mat.append_trait(mat_float64, trait = trait_object)
    assert numpy.all(mat.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(mat.trait == numpy.append(trait_object, trait_object, axis = 0))

### remove_trait

def test_remove_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "remove_trait")

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

### incorp_trait

def test_incorp_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "incorp_trait")

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

############################################################
##################### Sorting Methods ######################

### lexsort_trait

def test_lexsort_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "lexsort_trait")

def test_lexsort_trait_None(mat, trait_lexsort_indices):
    ix = mat.lexsort_trait(keys = None)
    assert numpy.all(ix == trait_lexsort_indices)

def test_lexsort_trait_tuple(mat, trait_object, trait_lexsort_indices):
    ix = mat.lexsort_trait(keys = (trait_object,))
    assert numpy.all(ix == trait_lexsort_indices)

### reorder_trait

def test_reorder_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "reorder_trait")

def test_reorder_trait_array_like(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.reorder_trait(trait_lexsort_indices)
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

### sort_trait

def test_sort_trait_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "sort_trait")

def test_sort_trait_None(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.sort_trait(keys = None)
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_sort_trait_tuple(mat, mat_float64, trait_object, trait_lexsort_indices):
    mat.sort_trait(keys = (trait_object,))
    assert numpy.all(mat.mat == mat_float64[trait_lexsort_indices])
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

############################################################
##################### Matrix File I/O ######################

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseTraitMatrix, "to_hdf5")

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

############################################################
##################### Matrix File I/O ######################

### from_hdf5

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseTraitMatrix, "from_hdf5")

def test_from_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    out = DenseTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_Path(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    out = DenseTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_h5py_File(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseTraitMatrix.from_hdf5(h5file)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    h5file.close()
    os.remove(fp)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseTraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTraitMatrix)

def test_check_is_DenseTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTraitMatrix(None, "mat")
