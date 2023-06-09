import pytest
import numpy
import copy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.DenseVariantMatrix import DenseVariantMatrix
from pybrops.core.mat.DenseVariantMatrix import check_is_DenseVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_int8():
    a = numpy.int8([[-1, 0, 1], [0, 1, -1], [1, -1, 0]])
    yield a

@pytest.fixture
def vrnt_chrgrp_int64():
    a = numpy.int64([1,2,2])
    yield a

@pytest.fixture
def vrnt_phypos_int64():
    a = numpy.int64([23, 7, 19])
    yield a

@pytest.fixture
def vrnt_name_object():
    a = numpy.object_(["snp1", "snp2", "snp3"])
    yield a

@pytest.fixture
def vrnt_genpos_float64():
    a = numpy.float64([0.3, 0.2, 0.4])
    yield a

@pytest.fixture
def vrnt_xoprob_float64():
    a = numpy.float64([0.5, 0.5, 0.3])
    yield a

@pytest.fixture
def vrnt_hapgrp_int64():
    a = numpy.int64([1,2,3])
    yield a

@pytest.fixture
def vrnt_hapalt_object():
    a = numpy.object_(["A", "T", "C"])
    yield a

@pytest.fixture
def vrnt_hapref_object():
    a = numpy.object_(["G", "A", "T"])
    yield a

@pytest.fixture
def vrnt_mask_bool():
    a = numpy.bool_([True, False, True])
    yield a

@pytest.fixture
def vrnt_chrgrp_name_int64():
    a = numpy.int64([1,2])
    yield a

@pytest.fixture
def vrnt_chrgrp_stix_int64():
    a = numpy.int64([0,1])
    yield a

@pytest.fixture
def vrnt_chrgrp_spix_int64():
    a = numpy.int64([1,3])
    yield a

@pytest.fixture
def vrnt_chrgrp_len_int64():
    a = numpy.int64([1,2])
    yield a

@pytest.fixture
def vrnt_lexsort_indices(vrnt_phypos_int64, vrnt_chrgrp_int64):
    a = numpy.lexsort((vrnt_phypos_int64, vrnt_chrgrp_int64))
    yield a

@pytest.fixture
def mat(mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    out = DenseVariantMatrix(
        mat = mat_int8,
        vrnt_chrgrp = vrnt_chrgrp_int64,
        vrnt_phypos = vrnt_phypos_int64,
        vrnt_name = vrnt_name_object,
        vrnt_genpos = vrnt_genpos_float64,
        vrnt_xoprob = vrnt_xoprob_float64,
        vrnt_hapgrp = vrnt_hapgrp_int64,
        vrnt_hapalt = vrnt_hapalt_object,
        vrnt_hapref = vrnt_hapref_object,
        vrnt_mask = vrnt_mask_bool
    )
    out.group_vrnt()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "__init__")

def test_copy_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "__copy__")

def test_deepcopy_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "__deepcopy__")

def test_adjoin_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "adjoin_vrnt")

def test_delete_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "delete_vrnt")

def test_insert_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "insert_vrnt")

def test_select_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "select_vrnt")

def test_concat_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "concat_vrnt")

def test_append_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "append_vrnt")

def test_remove_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "remove_vrnt")

def test_incorp_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "incorp_vrnt")

def test_lexsort_vrnt_is_concrete():
    assert_concrete_method(DenseVariantMatrix, "lexsort_vrnt")

# TODO: # FIXME: not_raises fails for an edge case
# def test_sort_vrnt_is_concrete():
#     generic_assert_abstract_method(DenseVariantMatrix, "sort_vrnt")
#
# def test_group_vrnt_is_concrete():
#     generic_assert_abstract_method(DenseVariantMatrix, "group_vrnt")
#
# def test_is_grouped_vrnt_is_concrete():
#     generic_assert_abstract_method(DenseVariantMatrix, "is_grouped_vrnt")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
def test_copy(mat):
    m = copy.copy(mat)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.vrnt_chrgrp == mat.vrnt_chrgrp)
    assert numpy.all(m.vrnt_phypos == mat.vrnt_phypos)
    assert numpy.all(m.vrnt_name == mat.vrnt_name)
    assert numpy.all(m.vrnt_genpos == mat.vrnt_genpos)
    assert numpy.all(m.vrnt_xoprob == mat.vrnt_xoprob)
    assert numpy.all(m.vrnt_hapgrp == mat.vrnt_hapgrp)
    assert numpy.all(m.vrnt_hapalt == mat.vrnt_hapalt)
    assert numpy.all(m.vrnt_hapref == mat.vrnt_hapref)
    assert numpy.all(m.vrnt_mask == mat.vrnt_mask)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_deepcopy(mat):
    m = copy.deepcopy(mat)
    # make sure object ID's are different
    assert id(m.mat) != id(mat.mat)
    assert id(m.vrnt_chrgrp) != id(mat.vrnt_chrgrp)
    assert id(m.vrnt_phypos) != id(mat.vrnt_phypos)
    assert id(m.vrnt_name) != id(mat.vrnt_name)
    assert id(m.vrnt_genpos) != id(mat.vrnt_genpos)
    assert id(m.vrnt_xoprob) != id(mat.vrnt_xoprob)
    assert id(m.vrnt_hapgrp) != id(mat.vrnt_hapgrp)
    assert id(m.vrnt_hapalt) != id(mat.vrnt_hapalt)
    assert id(m.vrnt_hapref) != id(mat.vrnt_hapref)
    assert id(m.vrnt_mask) != id(mat.vrnt_mask)
    assert id(m.vrnt_chrgrp_name) != id(mat.vrnt_chrgrp_name)
    assert id(m.vrnt_chrgrp_stix) != id(mat.vrnt_chrgrp_stix)
    assert id(m.vrnt_chrgrp_spix) != id(mat.vrnt_chrgrp_spix)
    assert id(m.vrnt_chrgrp_len) != id(mat.vrnt_chrgrp_len)
    # check that elements were copied correctly
    assert numpy.all(m.mat == mat.mat)
    assert numpy.all(m.vrnt_chrgrp == mat.vrnt_chrgrp)
    assert numpy.all(m.vrnt_phypos == mat.vrnt_phypos)
    assert numpy.all(m.vrnt_name == mat.vrnt_name)
    assert numpy.all(m.vrnt_genpos == mat.vrnt_genpos)
    assert numpy.all(m.vrnt_xoprob == mat.vrnt_xoprob)
    assert numpy.all(m.vrnt_hapgrp == mat.vrnt_hapgrp)
    assert numpy.all(m.vrnt_hapalt == mat.vrnt_hapalt)
    assert numpy.all(m.vrnt_hapref == mat.vrnt_hapref)
    assert numpy.all(m.vrnt_mask == mat.vrnt_mask)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

################################################################################
############################ Test Class Properties #############################
################################################################################

################### Taxa Data Properites ###################
def test_vrnt_chrgrp_fget(mat, vrnt_chrgrp_int64):
    assert numpy.all(mat.vrnt_chrgrp == vrnt_chrgrp_int64)

def test_vrnt_chrgrp_fset(mat, vrnt_chrgrp_int64):
    mat.vrnt_chrgrp = vrnt_chrgrp_int64
    assert numpy.all(mat.vrnt_chrgrp == vrnt_chrgrp_int64)

def test_vrnt_chrgrp_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_chrgrp

def test_vrnt_phypos_fget(mat, vrnt_phypos_int64):
    assert numpy.all(mat.vrnt_phypos == vrnt_phypos_int64)

def test_vrnt_phypos_fset(mat, vrnt_phypos_int64):
    mat.vrnt_phypos = vrnt_phypos_int64
    assert numpy.all(mat.vrnt_phypos == vrnt_phypos_int64)

def test_vrnt_phypos_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_phypos

def test_vrnt_name_fget(mat, vrnt_name_object):
    assert numpy.all(mat.vrnt_name == vrnt_name_object)

def test_vrnt_name_fset(mat, vrnt_name_object):
    mat.vrnt_name = vrnt_name_object
    assert numpy.all(mat.vrnt_name == vrnt_name_object)

def test_vrnt_name_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_name

def test_vrnt_genpos_fget(mat, vrnt_genpos_float64):
    assert numpy.all(mat.vrnt_genpos == vrnt_genpos_float64)

def test_vrnt_genpos_fset(mat, vrnt_genpos_float64):
    mat.vrnt_genpos = vrnt_genpos_float64
    assert numpy.all(mat.vrnt_genpos == vrnt_genpos_float64)

def test_vrnt_genpos_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_genpos

def test_vrnt_xoprob_fget(mat, vrnt_xoprob_float64):
    assert numpy.all(mat.vrnt_xoprob == vrnt_xoprob_float64)

def test_vrnt_xoprob_fset(mat, vrnt_xoprob_float64):
    mat.vrnt_xoprob = vrnt_xoprob_float64
    assert numpy.all(mat.vrnt_xoprob == vrnt_xoprob_float64)

def test_vrnt_xoprob_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_xoprob

def test_vrnt_hapgrp_fget(mat, vrnt_hapgrp_int64):
    assert numpy.all(mat.vrnt_hapgrp == vrnt_hapgrp_int64)

def test_vrnt_hapgrp_fset(mat, vrnt_hapgrp_int64):
    mat.vrnt_hapgrp = vrnt_hapgrp_int64
    assert numpy.all(mat.vrnt_hapgrp == vrnt_hapgrp_int64)

def test_vrnt_hapgrp_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_hapgrp

def test_vrnt_hapalt_fget(mat, vrnt_hapalt_object):
    assert numpy.all(mat.vrnt_hapalt == vrnt_hapalt_object)

def test_vrnt_hapalt_fset(mat, vrnt_hapalt_object):
    mat.vrnt_hapalt = vrnt_hapalt_object
    assert numpy.all(mat.vrnt_hapalt == vrnt_hapalt_object)

def test_vrnt_hapalt_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_hapalt

def test_vrnt_hapref_fget(mat, vrnt_hapref_object):
    assert numpy.all(mat.vrnt_hapref == vrnt_hapref_object)

def test_vrnt_hapref_fset(mat, vrnt_hapref_object):
    mat.vrnt_hapref = vrnt_hapref_object
    assert numpy.all(mat.vrnt_hapref == vrnt_hapref_object)

def test_vrnt_hapref_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_hapref

def test_vrnt_mask_fget(mat, vrnt_mask_bool):
    assert numpy.all(mat.vrnt_mask == vrnt_mask_bool)

def test_vrnt_mask_fset(mat, vrnt_mask_bool):
    mat.vrnt_mask = vrnt_mask_bool
    assert numpy.all(mat.vrnt_mask == vrnt_mask_bool)

def test_vrnt_mask_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_mask

################# Taxa Metadata Properites #################
def test_nvrnt_fget(mat, mat_int8):
    assert mat.nvrnt == len(mat_int8)

def test_nvrnt_fset(mat, mat_int8):
    with pytest.raises(AttributeError):
        mat.nvrnt = len(mat_int8)

def test_nvrnt_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.nvrnt

def test_vrnt_axis_fget(mat):
    assert mat.vrnt_axis == 0

def test_vrnt_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.vrnt_axis = 1

def test_vrnt_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_axis

def test_vrnt_chrgrp_name_fget(mat, vrnt_chrgrp_name_int64):
    assert numpy.all(mat.vrnt_chrgrp_name == vrnt_chrgrp_name_int64)

def test_vrnt_chrgrp_name_fset(mat, vrnt_chrgrp_name_int64):
    mat.vrnt_chrgrp_name = vrnt_chrgrp_name_int64
    assert numpy.all(mat.vrnt_chrgrp_name == vrnt_chrgrp_name_int64)

def test_vrnt_chrgrp_name_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_chrgrp_name

def test_vrnt_chrgrp_stix_fget(mat, vrnt_chrgrp_stix_int64):
    assert numpy.all(mat.vrnt_chrgrp_stix == vrnt_chrgrp_stix_int64)

def test_vrnt_chrgrp_stix_fset(mat, vrnt_chrgrp_stix_int64):
    mat.vrnt_chrgrp_stix = vrnt_chrgrp_stix_int64
    assert numpy.all(mat.vrnt_chrgrp_stix == vrnt_chrgrp_stix_int64)

def test_vrnt_chrgrp_stix_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_chrgrp_stix

def test_vrnt_chrgrp_spix_fget(mat, vrnt_chrgrp_spix_int64):
    assert numpy.all(mat.vrnt_chrgrp_spix == vrnt_chrgrp_spix_int64)

def test_vrnt_chrgrp_spix_fset(mat, vrnt_chrgrp_spix_int64):
    mat.vrnt_chrgrp_spix = vrnt_chrgrp_spix_int64
    assert numpy.all(mat.vrnt_chrgrp_spix == vrnt_chrgrp_spix_int64)

def test_vrnt_chrgrp_spix_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_chrgrp_spix

def test_vrnt_chrgrp_len_fget(mat, vrnt_chrgrp_len_int64):
    assert numpy.all(mat.vrnt_chrgrp_len == vrnt_chrgrp_len_int64)

def test_vrnt_chrgrp_len_fset(mat, vrnt_chrgrp_len_int64):
    mat.vrnt_chrgrp_len = vrnt_chrgrp_len_int64
    assert numpy.all(mat.vrnt_chrgrp_len == vrnt_chrgrp_len_int64)

def test_vrnt_chrgrp_len_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_chrgrp_len

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########### Matrix element copy-on-manipulation ############
def test_adjoin_vrnt_cls(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    m = mat.adjoin_vrnt(mat)
    assert numpy.all(m.mat == numpy.append(mat_int8, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.append(vrnt_chrgrp_int64, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.append(vrnt_phypos_int64, vrnt_phypos_int64, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.append(vrnt_name_object, vrnt_name_object, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.append(vrnt_genpos_float64, vrnt_genpos_float64, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.append(vrnt_xoprob_float64, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.append(vrnt_hapgrp_int64, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.append(vrnt_hapalt_object, vrnt_hapalt_object, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.append(vrnt_hapref_object, vrnt_hapref_object, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.append(vrnt_mask_bool, vrnt_mask_bool, axis = 0))

def test_adjoin_vrnt_ndarray(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    m = mat.adjoin_vrnt(
        mat_int8,
        vrnt_chrgrp = vrnt_chrgrp_int64,
        vrnt_phypos = vrnt_phypos_int64,
        vrnt_name = vrnt_name_object,
        vrnt_genpos = vrnt_genpos_float64,
        vrnt_xoprob = vrnt_xoprob_float64,
        vrnt_hapgrp = vrnt_hapgrp_int64,
        vrnt_hapalt = vrnt_hapalt_object,
        vrnt_hapref = vrnt_hapref_object,
        vrnt_mask = vrnt_mask_bool
    )
    assert numpy.all(m.mat == numpy.append(mat_int8, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.append(vrnt_chrgrp_int64, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.append(vrnt_phypos_int64, vrnt_phypos_int64, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.append(vrnt_name_object, vrnt_name_object, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.append(vrnt_genpos_float64, vrnt_genpos_float64, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.append(vrnt_xoprob_float64, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.append(vrnt_hapgrp_int64, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.append(vrnt_hapalt_object, vrnt_hapalt_object, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.append(vrnt_hapref_object, vrnt_hapref_object, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.append(vrnt_mask_bool, vrnt_mask_bool, axis = 0))

def test_delete_vrnt_cls_slice(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = slice(0,2,None)
    m = mat.delete_vrnt(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_delete_vrnt_cls_int(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = 0
    m = mat.delete_vrnt(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_delete_vrnt_cls_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [0,1,2]
    m = mat.delete_vrnt(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_insert_vrnt_cls_slice(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = slice(0,len(mat_int8),None)
    m = mat.insert_vrnt(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

def test_insert_vrnt_cls_int(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = 1
    m = mat.insert_vrnt(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

def test_insert_vrnt_cls_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [e for e in range(len(mat_int8))]
    m = mat.insert_vrnt(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

def test_select_vrnt_cls_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [0,0,1]
    m = mat.select_vrnt(obj)
    assert numpy.all(m.mat == numpy.take(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.take(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.take(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_name == numpy.take(vrnt_name_object, obj, axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.take(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.take(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.take(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.take(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.take(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.take(vrnt_mask_bool, obj, axis = 0))

def test_concat_vrnt_cls(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [mat, mat]
    m = mat.concat_vrnt(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_int8,mat_int8], axis = mat.vrnt_axis))
    assert numpy.all(m.vrnt_chrgrp == numpy.concatenate([vrnt_chrgrp_int64, vrnt_chrgrp_int64], axis = 0))
    assert numpy.all(m.vrnt_phypos == numpy.concatenate([vrnt_phypos_int64, vrnt_phypos_int64], axis = 0))
    assert numpy.all(m.vrnt_name == numpy.concatenate([vrnt_name_object, vrnt_name_object], axis = 0))
    assert numpy.all(m.vrnt_genpos == numpy.concatenate([vrnt_genpos_float64, vrnt_genpos_float64], axis = 0))
    assert numpy.all(m.vrnt_xoprob == numpy.concatenate([vrnt_xoprob_float64, vrnt_xoprob_float64], axis = 0))
    assert numpy.all(m.vrnt_hapgrp == numpy.concatenate([vrnt_hapgrp_int64, vrnt_hapgrp_int64], axis = 0))
    assert numpy.all(m.vrnt_hapalt == numpy.concatenate([vrnt_hapalt_object, vrnt_hapalt_object], axis = 0))
    assert numpy.all(m.vrnt_hapref == numpy.concatenate([vrnt_hapref_object, vrnt_hapref_object], axis = 0))
    assert numpy.all(m.vrnt_mask == numpy.concatenate([vrnt_mask_bool, vrnt_mask_bool], axis = 0))

########### Matrix element in-place-manipulation ###########
def test_append_vrnt_cls(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    mat.append_vrnt(mat)
    assert numpy.all(mat.mat == numpy.append(mat_int8, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.append(vrnt_chrgrp_int64, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.append(vrnt_phypos_int64, vrnt_phypos_int64, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.append(vrnt_name_object, vrnt_name_object, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.append(vrnt_genpos_float64, vrnt_genpos_float64, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.append(vrnt_xoprob_float64, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.append(vrnt_hapgrp_int64, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.append(vrnt_hapalt_object, vrnt_hapalt_object, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.append(vrnt_hapref_object, vrnt_hapref_object, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.append(vrnt_mask_bool, vrnt_mask_bool, axis = 0))

def test_append_vrnt_ndarray(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    mat.append_vrnt(
        mat_int8,
        vrnt_chrgrp = vrnt_chrgrp_int64,
        vrnt_phypos = vrnt_phypos_int64,
        vrnt_name = vrnt_name_object,
        vrnt_genpos = vrnt_genpos_float64,
        vrnt_xoprob = vrnt_xoprob_float64,
        vrnt_hapgrp = vrnt_hapgrp_int64,
        vrnt_hapalt = vrnt_hapalt_object,
        vrnt_hapref = vrnt_hapref_object,
        vrnt_mask = vrnt_mask_bool
    )
    assert numpy.all(mat.mat == numpy.append(mat_int8, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.append(vrnt_chrgrp_int64, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.append(vrnt_phypos_int64, vrnt_phypos_int64, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.append(vrnt_name_object, vrnt_name_object, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.append(vrnt_genpos_float64, vrnt_genpos_float64, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.append(vrnt_xoprob_float64, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.append(vrnt_hapgrp_int64, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.append(vrnt_hapalt_object, vrnt_hapalt_object, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.append(vrnt_hapref_object, vrnt_hapref_object, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.append(vrnt_mask_bool, vrnt_mask_bool, axis = 0))

def test_remove_vrnt_cls_slice(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = slice(0,2,None)
    mat.remove_vrnt(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_remove_vrnt_cls_int(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = 0
    mat.remove_vrnt(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_remove_vrnt_cls_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [0,1,2]
    mat.remove_vrnt(obj)
    assert numpy.all(mat.mat == numpy.delete(mat_int8, obj, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.delete(vrnt_chrgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.delete(vrnt_phypos_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.delete(vrnt_name_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.delete(vrnt_genpos_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.delete(vrnt_xoprob_float64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.delete(vrnt_hapgrp_int64, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.delete(vrnt_hapalt_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.delete(vrnt_hapref_object, obj, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.delete(vrnt_mask_bool, obj, axis = 0))

def test_incorp_vrnt_cls_slice(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = slice(0,len(mat_int8),None)
    mat.incorp_vrnt(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

def test_incorp_vrnt_cls_int(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = 1
    mat.incorp_vrnt(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

def test_incorp_vrnt_cls_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    obj = [e for e in range(len(mat_int8))]
    mat.incorp_vrnt(obj, mat)
    assert numpy.all(mat.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.vrnt_axis))
    assert numpy.all(mat.vrnt_chrgrp == numpy.insert(vrnt_chrgrp_int64, obj, vrnt_chrgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_phypos == numpy.insert(vrnt_phypos_int64, obj, vrnt_phypos_int64, axis = 0))
    assert numpy.all(mat.vrnt_name == numpy.insert(vrnt_name_object, obj, vrnt_name_object, axis = 0))
    assert numpy.all(mat.vrnt_genpos == numpy.insert(vrnt_genpos_float64, obj, vrnt_genpos_float64, axis = 0))
    assert numpy.all(mat.vrnt_xoprob == numpy.insert(vrnt_xoprob_float64, obj, vrnt_xoprob_float64, axis = 0))
    assert numpy.all(mat.vrnt_hapgrp == numpy.insert(vrnt_hapgrp_int64, obj, vrnt_hapgrp_int64, axis = 0))
    assert numpy.all(mat.vrnt_hapalt == numpy.insert(vrnt_hapalt_object, obj, vrnt_hapalt_object, axis = 0))
    assert numpy.all(mat.vrnt_hapref == numpy.insert(vrnt_hapref_object, obj, vrnt_hapref_object, axis = 0))
    assert numpy.all(mat.vrnt_mask == numpy.insert(vrnt_mask_bool, obj, vrnt_mask_bool, axis = 0))

##################### Sorting Methods ######################
def test_lexsort_vrnt_None(mat, vrnt_lexsort_indices):
    ix = mat.lexsort_vrnt(keys = None)
    assert numpy.all(ix == vrnt_lexsort_indices)

def test_lexsort_vrnt_tuple(mat, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_lexsort_indices):
    ix = mat.lexsort_vrnt(keys = (vrnt_phypos_int64, vrnt_chrgrp_int64))
    assert numpy.all(ix == vrnt_lexsort_indices)

def test_reorder_vrnt_array_like(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool, vrnt_lexsort_indices):
    mat.reorder_vrnt(vrnt_lexsort_indices)
    assert numpy.all(mat.mat == mat_int8[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_chrgrp == vrnt_chrgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_phypos == vrnt_phypos_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_name == vrnt_name_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_genpos == vrnt_genpos_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_xoprob == vrnt_xoprob_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapgrp == vrnt_hapgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapalt == vrnt_hapalt_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapref == vrnt_hapref_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_mask == vrnt_mask_bool[vrnt_lexsort_indices])

def test_sort_vrnt_None(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool, vrnt_lexsort_indices):
    mat.sort_vrnt(keys = None)
    assert numpy.all(mat.mat == mat_int8[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_chrgrp == vrnt_chrgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_phypos == vrnt_phypos_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_name == vrnt_name_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_genpos == vrnt_genpos_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_xoprob == vrnt_xoprob_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapgrp == vrnt_hapgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapalt == vrnt_hapalt_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapref == vrnt_hapref_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_mask == vrnt_mask_bool[vrnt_lexsort_indices])

def test_sort_vrnt_tuple(mat, mat_int8, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool, vrnt_lexsort_indices):
    mat.sort_vrnt(keys = (vrnt_phypos_int64, vrnt_chrgrp_int64))
    assert numpy.all(mat.mat == mat_int8[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_chrgrp == vrnt_chrgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_phypos == vrnt_phypos_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_name == vrnt_name_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_genpos == vrnt_genpos_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_xoprob == vrnt_xoprob_float64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapgrp == vrnt_hapgrp_int64[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapalt == vrnt_hapalt_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_hapref == vrnt_hapref_object[vrnt_lexsort_indices])
    assert numpy.all(mat.vrnt_mask == vrnt_mask_bool[vrnt_lexsort_indices])

def test_group_vrnt(mat, vrnt_chrgrp_name_int64, vrnt_chrgrp_stix_int64, vrnt_chrgrp_spix_int64, vrnt_chrgrp_len_int64):
    mat.group_vrnt()
    assert numpy.all(mat.vrnt_chrgrp_name == vrnt_chrgrp_name_int64)
    assert numpy.all(mat.vrnt_chrgrp_stix == vrnt_chrgrp_stix_int64)
    assert numpy.all(mat.vrnt_chrgrp_spix == vrnt_chrgrp_spix_int64)
    assert numpy.all(mat.vrnt_chrgrp_len == vrnt_chrgrp_len_int64)

def test_is_grouped_vrnt(mat):
    mat.group_vrnt()
    assert mat.is_grouped_vrnt() == (
        (mat.vrnt_chrgrp_name is not None) and
        (mat.vrnt_chrgrp_stix is not None) and
        (mat.vrnt_chrgrp_spix is not None) and
        (mat.vrnt_chrgrp_len is not None)
    )

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseVariantMatrix_is_concrete():
    assert_concrete_function(check_is_DenseVariantMatrix)

def test_check_is_DenseVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseVariantMatrix(None, "mat")
