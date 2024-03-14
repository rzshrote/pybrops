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

from pybrops.core.mat.DenseTaxaMatrix import DenseTaxaMatrix
from pybrops.core.mat.DenseTaxaMatrix import check_is_DenseTaxaMatrix

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
def taxa_grp_name_int64():
    a = numpy.int64([1,2])
    yield a

@pytest.fixture
def taxa_grp_stix_int64():
    a = numpy.int64([0,1])
    yield a

@pytest.fixture
def taxa_grp_spix_int64():
    a = numpy.int64([1,3])
    yield a

@pytest.fixture
def taxa_grp_len_int64():
    a = numpy.int64([1,2])
    yield a

@pytest.fixture
def taxa_lexsort_indices(taxa_object, taxa_grp_int64):
    a = numpy.lexsort((taxa_object, taxa_grp_int64))
    yield a

@pytest.fixture
def mat(mat_float64, taxa_object, taxa_grp_int64):
    out = DenseTaxaMatrix(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    out.group_taxa()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseTaxaMatrix)

################################################################################
######################## Test concrete special Methods #########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "__init__")

### __copy__

def test___copy___is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "__deepcopy__")


################################################################################
############################ Test Class Properties #############################
################################################################################

################### Taxa Data Properites ###################

### taxa

def test_taxa_fget(mat, taxa_object):
    assert numpy.all(mat.taxa == taxa_object)

def test_taxa_fset(mat, taxa_object):
    mat.taxa = taxa_object
    assert numpy.all(mat.taxa == taxa_object)

def test_taxa_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa

### taxa_grp

def test_taxa_grp_fget(mat, taxa_grp_int64):
    assert numpy.all(mat.taxa_grp == taxa_grp_int64)

def test_taxa_grp_fset(mat, taxa_grp_int64):
    mat.taxa_grp = taxa_grp_int64
    assert numpy.all(mat.taxa_grp == taxa_grp_int64)

def test_taxa_grp_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_grp

################# Taxa Metadata Properites #################

### ntaxa

def test_ntaxa_fget(mat, mat_float64):
    assert mat.ntaxa == len(mat_float64)

def test_ntaxa_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.ntaxa = len(mat_float64)

def test_ntaxa_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.ntaxa

### taxa_axis

def test_taxa_axis_fget(mat):
    assert mat.taxa_axis == 0

def test_taxa_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.taxa_axis = 1

def test_taxa_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_axis

### taxa_grp_name

def test_taxa_grp_name_fget(mat, taxa_grp_name_int64):
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)

def test_taxa_grp_name_fset(mat, taxa_grp_name_int64):
    mat.taxa_grp_name = taxa_grp_name_int64
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)

def test_taxa_grp_name_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_grp_name

### taxa_grp_stix

def test_taxa_grp_stix_fget(mat, taxa_grp_stix_int64):
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)

def test_taxa_grp_stix_fset(mat, taxa_grp_stix_int64):
    mat.taxa_grp_stix = taxa_grp_stix_int64
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)

def test_taxa_grp_stix_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_grp_stix

### taxa_grp_spix

def test_taxa_grp_spix_fget(mat, taxa_grp_spix_int64):
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)

def test_taxa_grp_spix_fset(mat, taxa_grp_spix_int64):
    mat.taxa_grp_spix = taxa_grp_spix_int64
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)

def test_taxa_grp_spix_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_grp_spix

### taxa_grp_len

def test_taxa_grp_len_fget(mat, taxa_grp_len_int64):
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)

def test_taxa_grp_len_fset(mat, taxa_grp_len_int64):
    mat.taxa_grp_len = taxa_grp_len_int64
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)

def test_taxa_grp_len_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_grp_len

################################################################################
############################# Test concrete methods ############################
################################################################################

### copy
        
def test_copy_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "copy")

def test_copy(mat):
    m = copy.copy(mat)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.taxa == mat.taxa)
    assert numpy.all(m.taxa_grp == mat.taxa_grp)
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

### deepcopy
        
def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "deepcopy")

def test_deepcopy(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = copy.deepcopy(mat)
    # make sure object ID's are different
    assert id(m.mat) != id(mat.mat)
    assert id(m.taxa) != id(mat.taxa)
    assert id(m.taxa_grp) != id(mat.taxa_grp)
    assert id(m.taxa_grp_name) != id(mat.taxa_grp_name)
    assert id(m.taxa_grp_stix) != id(mat.taxa_grp_stix)
    assert id(m.taxa_grp_spix) != id(mat.taxa_grp_spix)
    assert id(m.taxa_grp_len) != id(mat.taxa_grp_len)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.taxa == mat.taxa)
    assert numpy.all(m.taxa_grp == mat.taxa_grp)
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

############################################################
########### Matrix element copy-on-manipulation ############

### adjoin

def test_adjoin_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "adjoin")

### adjoin_taxa

def test_adjoin_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "adjoin_taxa")

def test_adjoin_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

def test_adjoin_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat_float64, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))

### delete

def test_delete_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "delete")

### delete_taxa

def test_delete_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "delete_taxa")

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

### insert

def test_insert_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "insert")

### insert_taxa

def test_insert_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "insert_taxa")

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

### select

def test_select_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "select")

### select_taxa

def test_select_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "select_taxa")

def test_select_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))

### concat

def test_concat_is_concrete():
    assert_classmethod_isconcrete(DenseTaxaMatrix, "concat")

### concat_taxa

def test_concat_taxa_is_concrete():
    assert_classmethod_isconcrete(DenseTaxaMatrix, "concat_taxa")

def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [mat, mat]
    m = mat.concat_taxa(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
    assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

############################################################
########### Matrix element in-place-manipulation ###########

### append

def test_append_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "append")

### append_taxa

def test_append_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "append_taxa")

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

### remove

def test_remove_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "remove")

### remove_taxa

def test_remove_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "remove_taxa")

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

### incorp

def test_incorp_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "incorp")

### incorp_taxa

def test_incorp_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "incorp_taxa")

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

############################################################
##################### Sorting Methods ######################

### lexsort

def test_lexsort_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "lexsort")

### lexsort_taxa

def test_lexsort_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "lexsort_taxa")

def test_lexsort_taxa_None(mat, taxa_lexsort_indices):
    ix = mat.lexsort_taxa(keys = None)
    assert numpy.all(ix == taxa_lexsort_indices)

def test_lexsort_taxa_tuple(mat, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    ix = mat.lexsort_taxa(keys = (taxa_object, taxa_grp_int64))
    assert numpy.all(ix == taxa_lexsort_indices)

### reorder

def test_reorder_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "reorder")

### reorder_taxa

def test_reorder_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "reorder_taxa")

def test_reorder_taxa_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices):
    mat.reorder_taxa(taxa_lexsort_indices)
    assert numpy.all(mat.mat == mat_float64[taxa_lexsort_indices])
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])

### sort

def test_sort_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "sort")

### sort_taxa

def test_sort_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "sort_taxa")

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
    assert_method_isconcrete(DenseTaxaMatrix, "group")

### group_taxa

def test_group_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "group_taxa")

def test_group_taxa(mat, taxa_grp_name_int64, taxa_grp_stix_int64, taxa_grp_spix_int64, taxa_grp_len_int64):
    mat.group_taxa()
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)
    assert mat.is_grouped_taxa()

### ungroup

def test_ungroup_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "ungroup")

### ungroup_taxa

def test_ungroup_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "ungroup_taxa")

def test_ungroup_taxa(mat):
    mat.ungroup_taxa()
    assert not mat.is_grouped_taxa()

### is_grouped

def test_is_grouped_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "is_grouped")

### is_grouped_taxa

def test_is_grouped_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaMatrix, "is_grouped_taxa")

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
    assert_method_isconcrete(DenseTaxaMatrix, "to_hdf5")

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
    assert_classmethod_isconcrete(DenseTaxaMatrix, "from_hdf5")

def test_from_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    out = DenseTaxaMatrix.from_hdf5(fp)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    os.remove(fp)

def test_from_hdf5_Path(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    out = DenseTaxaMatrix.from_hdf5(fp)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    os.remove(fp)

def test_from_hdf5_h5py_File(mat):
    fp = Path("tmp.h5")
    mat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseTaxaMatrix.from_hdf5(h5file)
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
    assert numpy.all(mat.taxa == out.taxa)
    assert numpy.all(mat.taxa_grp == out.taxa_grp)
    assert mat.ntaxa == out.ntaxa
    assert mat.taxa_axis == out.taxa_axis
    assert numpy.all(mat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(mat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(mat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(mat.taxa_grp_len == out.taxa_grp_len)
    h5file.close()
    os.remove(fp)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseTaxaMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTaxaMatrix)

def test_check_is_DenseTaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTaxaMatrix(None, "mat")
