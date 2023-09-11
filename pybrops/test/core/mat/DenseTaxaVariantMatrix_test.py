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

from pybrops.core.mat.DenseTaxaVariantMatrix import DenseTaxaVariantMatrix
from pybrops.core.mat.DenseTaxaVariantMatrix import check_is_DenseTaxaVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_int8():
    a = numpy.float64([[0, 1, 1], [0, 0, 1], [1, 1, 1]])
    yield a

###################### Taxa fixtures #######################
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

##################### Variant fixtures #####################
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

############################################################
@pytest.fixture
def mat(mat_int8, taxa_object, taxa_grp_int64, vrnt_chrgrp_int64, vrnt_phypos_int64, vrnt_name_object, vrnt_genpos_float64, vrnt_xoprob_float64, vrnt_hapgrp_int64, vrnt_hapalt_object, vrnt_hapref_object, vrnt_mask_bool):
    a = DenseTaxaVariantMatrix(
        mat = mat_int8,
        taxa = taxa_object,
        taxa_grp = taxa_grp_int64,
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
    a.group_taxa()
    a.group_vrnt()
    yield a

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseTaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "__init__")

def test_copy_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "__copy__")

def test_deepcopy_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "__deepcopy__")

def test_adjoin_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "adjoin_taxa")

def test_adjoin_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "adjoin_vrnt")

def test_delete_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "delete_taxa")

def test_delete_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "delete_vrnt")

def test_insert_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "insert_taxa")

def test_insert_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "insert_vrnt")

def test_select_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "select_taxa")

def test_select_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "select_vrnt")

def test_concat_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "concat_taxa")

def test_concat_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "concat_vrnt")

def test_append_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "append_taxa")

def test_append_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "append_vrnt")

def test_remove_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "remove_taxa")

def test_remove_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "remove_vrnt")

def test_incorp_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "incorp_taxa")

def test_incorp_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "incorp_vrnt")

def test_lexsort_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "lexsort_taxa")

def test_lexsort_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "lexsort_vrnt")

def test_sort_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "sort_taxa")

def test_sort_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "sort_vrnt")

def test_group_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "group_taxa")

def test_group_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "group_vrnt")

def test_ungroup_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "ungroup_taxa")

def test_ungroup_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "ungroup_vrnt")

def test_is_grouped_taxa_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "is_grouped_taxa")

def test_is_grouped_vrnt_is_concrete():
    assert_concrete_method(DenseTaxaVariantMatrix, "is_grouped_vrnt")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################
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

def test_deepcopy(mat, mat_int8, taxa_object, taxa_grp_int64):
    m = copy.deepcopy(mat)
    # make sure object ID's are different
    assert id(m.mat) != id(mat.mat)
    assert id(m.taxa) != id(mat.taxa)
    assert id(m.taxa_grp) != id(mat.taxa_grp)
    assert id(m.taxa_grp_name) != id(mat.taxa_grp_name)
    assert id(m.taxa_grp_stix) != id(mat.taxa_grp_stix)
    assert id(m.taxa_grp_spix) != id(mat.taxa_grp_spix)
    assert id(m.taxa_grp_len) != id(mat.taxa_grp_len)
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
    assert numpy.all(m.taxa == mat.taxa)
    assert numpy.all(m.taxa_grp == mat.taxa_grp)
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
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

################# Taxa Metadata Properites #################
def test_taxa_axis_fget(mat):
    assert mat.taxa_axis == 0

def test_taxa_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.taxa_axis = 1

def test_taxa_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_axis

############### Variant Metadata Properites ################
def test_vrnt_axis_fget(mat):
    assert mat.vrnt_axis == 1

def test_vrnt_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.vrnt_axis = 1

def test_vrnt_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_axis

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########### Matrix element copy-on-manipulation ############
def test_adjoin_taxa_cls(mat, mat_int8, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat)
    assert numpy.all(m.mat == numpy.append(mat_int8, mat_int8, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_adjoin_taxa_ndarray(mat, mat_int8, taxa_object, taxa_grp_int64):
    m = mat.adjoin_taxa(mat_int8, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    assert numpy.all(m.mat == numpy.append(mat_int8, mat_int8, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_delete_taxa_cls_slice(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = slice(0,2,None)
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_delete_taxa_cls_int(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = 0
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_delete_taxa_cls_array_like(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = [0,1,2]
    m = mat.delete_taxa(obj)
    assert numpy.all(m.mat == numpy.delete(mat_int8, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_insert_taxa_cls_slice(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = slice(0,len(mat_int8),None)
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_insert_taxa_cls_int(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = 1
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_insert_taxa_cls_array_like(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = [e for e in range(len(mat_int8))]
    m = mat.insert_taxa(obj, mat)
    assert numpy.all(m.mat == numpy.insert(mat_int8, obj, mat_int8, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.insert(taxa_object, obj, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.insert(taxa_grp_int64, obj, taxa_grp_int64, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_select_taxa_cls_array_like(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    assert numpy.all(m.mat == numpy.take(mat_int8, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_concat_taxa_cls(mat, mat_int8, taxa_object, taxa_grp_int64):
    obj = [mat, mat]
    m = mat.concat_taxa(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_int8,mat_int8], axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
    assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

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
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseTaxaVariantMatrix_is_concrete():
    assert_concrete_function(check_is_DenseTaxaVariantMatrix)

def test_check_is_DenseTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTaxaVariantMatrix(None, "mat")
