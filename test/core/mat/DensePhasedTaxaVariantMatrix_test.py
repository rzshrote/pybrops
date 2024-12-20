import pytest
import numpy

from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DensePhasedTaxaVariantMatrix import DensePhasedTaxaVariantMatrix
from pybrops.core.mat.DensePhasedTaxaVariantMatrix import check_is_DensePhasedTaxaVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_int8():
    a = numpy.int8([
        [[0, 1, 1], [0, 0, 1], [1, 1, 1]],
        [[1, 0, 0], [1, 1, 1], [0, 0, 0]],
        [[1, 0, 1], [1, 1, 0], [1, 1, 0]]
    ])
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
    yield DensePhasedTaxaVariantMatrix(
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

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DensePhasedTaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "__init__")

def test_adjoin_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "adjoin_phase")

def test_delete_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "delete_phase")

def test_insert_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "insert_phase")

def test_select_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "select_phase")

def test_concat_phase_is_concrete():
    assert_classmethod_isconcrete(DensePhasedTaxaVariantMatrix, "concat_phase")

def test_append_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "append_phase")

def test_remove_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "remove_phase")

def test_incorp_phase_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "incorp_phase")

def test_group_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "group")

def test_ungroup_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "ungroup")

def test_is_grouped_is_concrete():
    assert_method_isconcrete(DensePhasedTaxaVariantMatrix, "is_grouped")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################ Phase Metadata Properites #################
def test_phase_axis_fget(mat):
    assert mat.phase_axis == 0

def test_phase_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.phase_axis = 1

def test_phase_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.phase_axis

################# Taxa Metadata Properites #################
def test_taxa_axis_fget(mat):
    assert mat.taxa_axis == 1

def test_taxa_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.taxa_axis = 1

def test_taxa_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.taxa_axis

############### Variant Metadata Properites ################
def test_vrnt_axis_fget(mat):
    assert mat.vrnt_axis == 2

def test_vrnt_axis_fset(mat):
    with pytest.raises(AttributeError):
        mat.vrnt_axis = 1

def test_vrnt_axis_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.vrnt_axis

################################################################################
###################### Test concrete method functionality ######################
################################################################################

### adjoin ###
def test_adjoin_phase_cls(mat, mat_int8):
    m = mat.adjoin_phase(mat)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_adjoin_phase_ndarray(mat, mat_int8):
    m = mat.adjoin_phase(mat_int8)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

### adjoin_phase ###

### delete ###
def test_delete_phase_cls_slice(mat, mat_int8):
    obj = slice(0,2,None)
    m = mat.delete_phase(obj)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_delete_phase_cls_int(mat, mat_int8):
    obj = 0
    m = mat.delete_phase(obj)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_delete_phase_cls_array_like(mat, mat_int8):
    obj = [0,1,2]
    m = mat.delete_phase(obj)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

### delete_phase ###

### insert ###
def test_insert_phase_cls_slice(mat, mat_int8):
    obj = slice(0,len(mat_int8),None)
    m = mat.insert_phase(obj, mat)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_insert_phase_cls_int(mat, mat_int8):
    obj = 1
    m = mat.insert_phase(obj, mat)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

def test_insert_phase_cls_array_like(mat, mat_int8):
    obj = [e for e in range(len(mat_int8))]
    m = mat.insert_phase(obj, mat)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

### insert_phase ###

### select ###
def test_select_phase_cls_array_like(mat, mat_int8):
    obj = [0,0,1]
    m = mat.select_phase(obj)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

### select_phase ###

### concat ###
def test_concat_phase_cls(mat, mat_int8):
    obj = [mat, mat]
    m = mat.concat_phase(obj)
    # test metadata copying
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)
    assert numpy.all(m.vrnt_chrgrp_name == mat.vrnt_chrgrp_name)
    assert numpy.all(m.vrnt_chrgrp_stix == mat.vrnt_chrgrp_stix)
    assert numpy.all(m.vrnt_chrgrp_spix == mat.vrnt_chrgrp_spix)
    assert numpy.all(m.vrnt_chrgrp_len == mat.vrnt_chrgrp_len)

### concat_phase ###

### append ###

### remove ###

### incorp ###

### group ###
def test_group_axis_0(mat):
    with pytest.raises(ValueError):
        mat.ungroup(0)
    with pytest.raises(ValueError):
        mat.group(0)
    assert not mat.is_grouped(0)

def test_group_axis_1(mat):
    with not_raises(ValueError):
        mat.ungroup(1)
    with not_raises(ValueError):
        mat.group(1)
    assert mat.is_grouped(1)

def test_group_axis_2(mat):
    with not_raises(ValueError):
        mat.ungroup(2)
    with not_raises(ValueError):
        mat.group(2)
    assert mat.is_grouped(2)

### ungroup ###
def test_ungroup_axis_0(mat):
    with pytest.raises(ValueError):
        mat.ungroup(0)
    assert not mat.is_grouped(0)
    with pytest.raises(ValueError):
        mat.group(0)
    assert not mat.is_grouped(0)

def test_ungroup_axis_1(mat):
    with not_raises(ValueError):
        mat.ungroup(1)
    assert not mat.is_grouped(1)
    with not_raises(ValueError):
        mat.group(1)
    assert mat.is_grouped(1)

def test_ungroup_axis_2(mat):
    with not_raises(ValueError):
        mat.ungroup(2)
    assert not mat.is_grouped(2)
    with not_raises(ValueError):
        mat.group(2)
    assert mat.is_grouped(2)

### is_grouped ###
def test_is_grouped_axis_0(mat):
    with pytest.raises(ValueError):
        mat.ungroup(0)
    value1 = mat.is_grouped(0)
    assert not value1
    with pytest.raises(ValueError):
        mat.group(0)
    value2 = mat.is_grouped(0)
    assert not value2
    assert value1 == False
    assert value2 == False
    assert value1 == value2

def test_is_grouped_axis_1(mat):
    with not_raises(ValueError):
        mat.ungroup(1)
    value1 = mat.is_grouped(1)
    assert not value1
    with not_raises(ValueError):
        mat.group(1)
    value2 = mat.is_grouped(1)
    assert value2
    assert value1 == False
    assert value2 == True
    assert value1 != value2

def test_is_grouped_axis_2(mat):
    with not_raises(ValueError):
        mat.ungroup(2)
    value1 = mat.is_grouped(2)
    assert not value1
    with not_raises(ValueError):
        mat.group(2)
    value2 = mat.is_grouped(2)
    assert value2
    assert value1 == False
    assert value2 == True
    assert value1 != value2

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DensePhasedTaxaVariantMatrix_is_concrete():
    assert_function_isconcrete(check_is_DensePhasedTaxaVariantMatrix)

def test_check_is_DensePhasedTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_DensePhasedTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DensePhasedTaxaVariantMatrix(None, "mat")
