import pytest
import numpy
import copy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import DenseTaxaTraitMatrix
from pybropt.core.mat import is_DenseTaxaTraitMatrix
from pybropt.core.mat import check_is_DenseTaxaTraitMatrix
from pybropt.core.mat import cond_check_is_DenseTaxaTraitMatrix

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
    generic_assert_docstring(DenseTaxaTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "__init__")

def test_copy_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "__copy__")

def test_deepcopy_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "__deepcopy__")

def test_adjoin_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "adjoin_taxa")

def test_adjoin_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "adjoin_trait")

def test_delete_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "delete_taxa")

def test_delete_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "delete_trait")

def test_insert_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "insert_taxa")

def test_insert_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "insert_trait")

def test_select_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "select_taxa")

def test_select_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "select_trait")

def test_concat_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "concat_taxa")

def test_concat_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "concat_trait")

def test_append_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "append_taxa")

def test_append_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "append_trait")

def test_remove_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "remove_taxa")

def test_remove_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "remove_trait")

def test_incorp_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "incorp_taxa")

def test_incorp_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "incorp_trait")

def test_lexsort_taxa_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "lexsort_taxa")

def test_lexsort_trait_is_concrete():
    generic_assert_concrete_method(DenseTaxaTraitMatrix, "lexsort_trait")

# TODO: # FIXME: not_raises fails for an edge case
# def test_sort_taxa_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "sort_taxa")
#
# def test_sort_trait_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "sort_trait")
#
# def test_group_taxa_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "group_taxa")
#
# def test_group_trait_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "group_trait")
#
# def test_is_grouped_taxa_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "is_grouped_taxa")
#
# def test_is_grouped_trait_is_concrete():
#     generic_assert_abstract_method(DenseTaxaTraitMatrix, "is_grouped_trait")

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
    assert numpy.all(m.trait == mat.trait)

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

################# Taxa Metadata Properites #################
def test_trait_axis_fget(mat):
    assert mat.trait_axis == 1

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

def test_select_taxa_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [0,0,1]
    m = mat.select_taxa(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))

def test_select_trait_cls_array_like(mat, mat_float64, trait_object):
    obj = [0,0,1]
    m = mat.select_trait(obj)
    assert numpy.all(m.mat == numpy.take(mat_float64, obj, axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
    obj = [mat, mat]
    m = mat.concat_taxa(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
    assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
    assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

def test_concat_trait_cls(mat, mat_float64, trait_object):
    obj = [mat, mat]
    m = mat.concat_trait(obj)
    assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.trait_axis))
    assert numpy.all(m.trait == numpy.concatenate([trait_object, trait_object], axis = 0))
    assert numpy.all(m.taxa_grp_name == mat.taxa_grp_name)
    assert numpy.all(m.taxa_grp_stix == mat.taxa_grp_stix)
    assert numpy.all(m.taxa_grp_spix == mat.taxa_grp_spix)
    assert numpy.all(m.taxa_grp_len == mat.taxa_grp_len)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseTaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseTaxaTraitMatrix)

def test_is_DenseTaxaTraitMatrix(mat):
    assert is_DenseTaxaTraitMatrix(mat)

def test_check_is_DenseTaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseTaxaTraitMatrix)

def test_check_is_DenseTaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseTaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseTaxaTraitMatrix(None, "mat")

def test_cond_check_is_DenseTaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseTaxaTraitMatrix)
