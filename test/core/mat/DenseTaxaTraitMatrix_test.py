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

from pybrops.core.mat.DenseTaxaTraitMatrix import DenseTaxaTraitMatrix
from pybrops.core.mat.DenseTaxaTraitMatrix import check_is_DenseTaxaTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([[3.3, 9.2, 5.6], [8.7, 3.7, 4.1], [9.0, 4.7, 3.8]])
    yield a

### Taxa fixtures
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

### Trait fixtures
@pytest.fixture
def trait_object():
    a = numpy.object_(["yield", "oil", "protein"])
    yield a

@pytest.fixture
def trait_lexsort_indices(trait_object):
    a = numpy.lexsort((trait_object,))
    yield a

@pytest.fixture
def mat(mat_float64, taxa_object, taxa_grp_int64, trait_object):
    out = DenseTaxaTraitMatrix(
        mat = mat_float64,
        taxa = taxa_object,
        taxa_grp = taxa_grp_int64,
        trait = trait_object
    )
    out.group_taxa()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseTaxaTraitMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "__init__")

### __copy__

def test___copy___is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "__deepcopy__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################# Taxa Metadata Properites #################

### taxa_axis

def test_taxa_axis_fget(mat):
    assert mat.taxa_axis == 0

def test_taxa_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.taxa_axis = 1

def test_taxa_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_axis

################# Taxa Metadata Properites #################

### trait_axis

def test_trait_axis_fget(mat):
    assert mat.trait_axis == 1

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
    assert numpy.all(m.trait == mat.trait)

### deepcopy

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
    assert id(m.trait) != id(mat.trait)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.taxa == mat.taxa)
    assert numpy.all(m.taxa_grp == mat.taxa_grp)
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.trait == mat.trait)

############################################################
########### Matrix element copy-on-manipulation ############

### adjoin

def test_adjoin_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "adjoin")

### adjoin_taxa

def test_adjoin_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "adjoin_taxa")

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

### adjoin_trait

def test_adjoin_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "adjoin_trait")

def test_adjoin_trait_cls(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_adjoin_trait_ndarray(mat, mat_float64, trait_object):
    m = mat.adjoin_trait(mat_float64, trait = trait_object)
    assert numpy.all(m.mat == numpy.append(mat_float64, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

### delete

def test_delete_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "delete")

### delete_taxa

def test_delete_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "delete_taxa")

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

### delete_trait

def test_delete_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "delete_trait")

def test_delete_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,2,None)
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_delete_trait_cls_int(mat, mat_float64, trait_object):
    obj = 0
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_delete_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,1,2]
    m = mat.delete_trait(obj)
    assert numpy.all(m.mat == numpy.delete(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

### insert

def test_insert_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "insert")

### insert_taxa

def test_insert_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "insert_taxa")

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

### insert_trait

def test_insert_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "insert_trait")

def test_insert_trait_cls_slice(mat, mat_float64, trait_object):
    obj = slice(0,len(mat_float64),None)
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_insert_trait_cls_int(mat, mat_float64, trait_object):
    obj = 1
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_insert_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [e for e in range(len(mat_float64))]
    m = mat.insert_trait(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_float64, obj, mat_float64, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.insert(trait_object, obj, trait_object, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

### select

def test_select_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "select")

### select_taxa

def test_select_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "select_taxa")

def test_select_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))

### select_trait

def test_select_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "select_trait")

def test_select_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,0,1]
    m = mat.select_trait(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

### concat

def test_concat_is_concrete():
    assert_classmethod_isconcrete(DenseTaxaTraitMatrix, "concat")

### concat_taxa

def test_concat_taxa_is_concrete():
    assert_classmethod_isconcrete(DenseTaxaTraitMatrix, "concat_taxa")

def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [mat, mat]
    m = mat.concat_taxa(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
    assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

### concat_trait

def test_concat_trait_is_concrete():
    assert_classmethod_isconcrete(DenseTaxaTraitMatrix, "concat_trait")

def test_concat_trait_cls(mat, mat_float64, trait_object):
    obj = [mat, mat]
    m = mat.concat_trait(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.concatenate([trait_object, trait_object], axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

############################################################
########### Matrix element in-place-manipulation ###########

### append

def test_append_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "append")

### append_taxa

def test_append_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "append_taxa")

### append_trait

def test_append_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "append_trait")

### remove

def test_remove_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "remove")

### remove_taxa

def test_remove_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "remove_taxa")

### remove_trait

def test_remove_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "remove_trait")

### incorp

def test_incorp_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "incorp")

### incorp_taxa

def test_incorp_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "incorp_taxa")

### incorp_trait

def test_incorp_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "incorp_trait")

############################################################
##################### Sorting Methods ######################

### lexsort

def test_lexsort_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "lexsort")

### lexsort_taxa

def test_lexsort_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "lexsort_taxa")

### lexsort_trait

def test_lexsort_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "lexsort_trait")

### sort

def test_sort_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "sort")

### sort_taxa

def test_sort_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "sort_taxa")

### sort_trait

def test_sort_trait_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "sort_trait")

############################################################
##################### Grouping Methods #####################

### group

def test_group_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "group")

### group_taxa

def test_group_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "group_taxa")

### is_grouped

def test_is_grouped_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "is_grouped")

### is_grouped_taxa

def test_is_grouped_taxa_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "is_grouped_taxa")

############################################################
##################### Matrix File I/O ######################

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseTaxaTraitMatrix, "to_hdf5")

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
    assert_classmethod_isconcrete(DenseTaxaTraitMatrix, "from_hdf5")

def test_from_hdf5_str(mat):
    fp = "tmp.h5"
    mat.to_hdf5(fp)
    out = DenseTaxaTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
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
    out = DenseTaxaTraitMatrix.from_hdf5(fp)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
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
    out = DenseTaxaTraitMatrix.from_hdf5(h5file)
    # general
    assert numpy.all(mat.mat == out.mat)
    assert mat.mat_ndim == out.mat_ndim
    assert mat.mat_shape == out.mat_shape
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
def test_check_is_DenseTaxaTraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTaxaTraitMatrix)

def test_check_is_DenseTaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTaxaTraitMatrix(None, "mat")
