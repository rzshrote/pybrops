import os
from pathlib import Path
import pytest
import numpy
import h5py

from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import check_is_DenseSquareTaxaTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [[3.3, 9.2, 5.6, 2.5],
         [8.7, 3.7, 4.1, 6.2],
         [9.0, 4.7, 3.8, 1.5],
         [0.6, 5.2, 8.3, 1.8]],
        [[2.2, 3.2, 5.2, 7.2],
         [2.3, 3.3, 5.3, 7.3],
         [2.5, 3.5, 5.5, 7.5],
         [2.7, 3.7, 5.7, 7.7]], 
    ]).transpose(1,2,0)
    yield a

@pytest.fixture
def taxa_object():
    a = numpy.object_(["A", "B", "C", "D"])
    yield a

@pytest.fixture
def taxa_grp_int64():
    a = numpy.int64([1,1,2,2])
    yield a

@pytest.fixture
def taxa_grp_name_int64():
    a = numpy.int64([1,2])
    yield a

@pytest.fixture
def taxa_grp_stix_int64():
    a = numpy.int64([0,2])
    yield a

@pytest.fixture
def taxa_grp_spix_int64():
    a = numpy.int64([2,4])
    yield a

@pytest.fixture
def taxa_grp_len_int64():
    a = numpy.int64([2,2])
    yield a

@pytest.fixture
def taxa_lexsort_indices(taxa_object, taxa_grp_int64):
    a = numpy.lexsort((taxa_object, taxa_grp_int64))
    yield a

@pytest.fixture
def mat(mat_float64, taxa_object, taxa_grp_int64):
    out = DenseSquareTaxaTraitMatrix(
        mat = mat_float64,
        taxa = taxa_object,
        taxa_grp = taxa_grp_int64
    )
    out.group_taxa()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseSquareTaxaTraitMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "__init__")

def test___copy___is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "__copy__")

def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "__deepcopy__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_copy_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "copy")

def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "deepcopy")

############################################################
########### Matrix element copy-on-manipulation ############

### adjoin

def test_adjoin_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "adjoin")

def test_adjoin_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    upper = numpy.concatenate([a,b], axis = 1)
    lower = numpy.concatenate([b,a], axis = 1)
    mtrue = numpy.concatenate([upper,lower], axis = 0)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_adjoin_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    upper = numpy.concatenate([a,b], axis = 1)
    lower = numpy.concatenate([b,a], axis = 1)
    mtrue = numpy.concatenate([upper,lower], axis = 0)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

### delete

def test_delete_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "delete")

def test_delete_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,2,None)
    m = mat.delete_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_delete_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 0
    m = mat.delete_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_delete_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,1,2]
    m = mat.delete_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

### insert

def test_insert_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "insert")

# def test_insert_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = slice(0,len(mat_float64),None)
#     m = mat.insert_taxa(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
#
# def test_insert_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = 1
#     m = mat.insert_taxa(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
#
# def test_insert_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = [e for e in range(len(mat_float64))]
#     m = mat.insert_taxa(obj, mat)
#     assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

### select

def test_select_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "select")

def test_select_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.take(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))

### concat

def test_concat_is_concrete():
    assert_classmethod_isconcrete(DenseSquareTaxaTraitMatrix, "concat")

# def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = [mat, mat]
#     m = mat.concat_taxa(obj)
#     assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
#     assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
#     assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

########### Matrix element in-place-manipulation ###########

### append

def test_append_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "append")

def test_append_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    mat.append_taxa(mat)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    upper = numpy.concatenate([a,b], axis = 1)
    lower = numpy.concatenate([b,a], axis = 1)
    mtrue = numpy.concatenate([upper,lower], axis = 0)
    assert numpy.array_equal(mat.mat, mtrue, equal_nan = True)
    assert numpy.all(mat.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_append_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64):
    mat.append_taxa(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    upper = numpy.concatenate([a,b], axis = 1)
    lower = numpy.concatenate([b,a], axis = 1)
    mtrue = numpy.concatenate([upper,lower], axis = 0)
    assert numpy.array_equal(mat.mat, mtrue, equal_nan = True)
    assert numpy.all(mat.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

### remove

def test_remove_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "remove")

def test_remove_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = slice(0,2,None)
    mat.remove_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_remove_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = 0
    mat.remove_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

def test_remove_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,1,2]
    mat.remove_taxa(obj)
    mtrue = mat_float64.copy()
    for axis in mat.square_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(mat.mat == mtrue)
    assert numpy.all(mat.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(mat.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))

### incorp

def test_incorp_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "incorp")

# def test_incorp_taxa_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = slice(0,len(mat_float64),None)
#     mat.incorp_taxa(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
#
# def test_incorp_taxa_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = 1
#     mat.incorp_taxa(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
#
# def test_incorp_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = [e for e in range(len(mat_float64))]
#     mat.incorp_taxa(obj, mat)
#     assert numpy.all(mat.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.taxa_axis))
#     assert numpy.all(mat.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
#     assert numpy.all(mat.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))

##################### Sorting Methods ######################

### lexsort

def test_lexsort_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "lexsort")

def test_lexsort_taxa_None(mat, taxa_lexsort_indices):
    ix = mat.lexsort_taxa(keys = None)
    assert numpy.all(ix == taxa_lexsort_indices)

def test_lexsort_taxa_tuple(mat, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    ix = mat.lexsort_taxa(keys = (taxa_object, taxa_grp_int64))
    assert numpy.all(ix == taxa_lexsort_indices)

### reorder

def test_reorder_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "reorder")

def test_reorder_taxa_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    mat.reorder_taxa(taxa_lexsort_indices)
    assert numpy.all(mat.mat == mat_float64[taxa_lexsort_indices])
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])

### sort

def test_sort_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "sort")

def test_sort_taxa_None(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    mat.sort_taxa(keys = None)
    assert numpy.all(mat.mat == mat_float64[taxa_lexsort_indices])
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])

def test_sort_taxa_tuple(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    mat.sort_taxa(keys = (taxa_object, taxa_grp_int64))
    assert numpy.all(mat.mat == mat_float64[taxa_lexsort_indices])
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])

### group

def test_group_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "group")

def test_group_taxa(mat, taxa_grp_name_int64, taxa_grp_stix_int64, taxa_grp_spix_int64, taxa_grp_len_int64):
    mat.group_taxa()
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)

### ungroup

def test_ungroup_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "ungroup")

### is_grouped

def test_is_grouped_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "is_grouped")

def test_is_grouped_taxa(mat):
    mat.group_taxa()
    assert mat.is_grouped_taxa() == (
        (mat.taxa_grp_name is not None) and
        (mat.taxa_grp_stix is not None) and
        (mat.taxa_grp_spix is not None) and
        (mat.taxa_grp_len is not None)
    )

############################################################
##################### Matrix File I/O ######################

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseSquareTaxaTraitMatrix, "to_hdf5")

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
    assert_classmethod_isconcrete(DenseSquareTaxaTraitMatrix, "from_hdf5")

def test_from_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    out = DenseSquareTaxaTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert mat.nsquare == out.nsquare
    assert mat.square_axes == out.square_axes
    assert mat.square_axes_len == out.square_axes_len
    assert mat.square_taxa_axes == out.square_taxa_axes
    # taxa
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_Path(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    out = DenseSquareTaxaTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert mat.nsquare == out.nsquare
    assert mat.square_axes == out.square_axes
    assert mat.square_axes_len == out.square_axes_len
    assert mat.square_taxa_axes == out.square_taxa_axes
    # taxa
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_h5py_File(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseSquareTaxaTraitMatrix.from_hdf5(h5file)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert mat.nsquare == out.nsquare
    assert mat.square_axes == out.square_axes
    assert mat.square_axes_len == out.square_axes_len
    assert mat.square_taxa_axes == out.square_taxa_axes
    # taxa
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(mat.trait == out.trait)
    assert mat.ntrait == out.ntrait
    assert mat.trait_axis == out.trait_axis
    h5file.close()
    os.remove(fp)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSquareTaxaTraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseSquareTaxaTraitMatrix)

def test_check_is_DenseSquareTaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseSquareTaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseSquareTaxaTraitMatrix(None, "mat")
