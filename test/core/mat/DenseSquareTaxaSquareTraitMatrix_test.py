import pytest
import numpy
import copy

from pybrops.test.assert_python import assert_concrete_property, not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.DenseSquareTaxaSquareTraitMatrix import DenseSquareTaxaSquareTraitMatrix
from pybrops.core.mat.DenseSquareTaxaSquareTraitMatrix import check_is_DenseSquareTaxaSquareTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ntaxa():
    yield 6 # must be divisible by 2

@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def mat_float64(ntaxa,ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntrait,ntrait))
    yield out

@pytest.fixture
def taxa_object(ntaxa):
    out = numpy.array(["Taxon"+str(i).zfill(3) for i in range(ntaxa)], dtype = object)
    yield out

@pytest.fixture
def trait_object(ntrait):
    out = numpy.array(["Trait"+str(i) for i in range(ntrait)], dtype = object)
    yield out

@pytest.fixture
def taxa_grp_int64(ntaxa):
    out = numpy.repeat(numpy.arange(ntaxa//2)+1, 2)
    yield out

@pytest.fixture
def taxa_grp_name_int64(ntaxa):
    out = numpy.arange(ntaxa//2)+1
    yield out

@pytest.fixture
def taxa_grp_stix_int64(ntaxa):
    out = numpy.arange(0, ntaxa, 2)
    yield out

@pytest.fixture
def taxa_grp_spix_int64(ntaxa):
    out = numpy.arange(0, ntaxa, 2) + 2
    yield out

@pytest.fixture
def taxa_grp_len_int64(ntaxa):
    out = numpy.repeat(2, ntaxa//2)
    yield out

@pytest.fixture
def taxa_lexsort_indices(taxa_object, taxa_grp_int64):
    a = numpy.lexsort((taxa_object, taxa_grp_int64))
    yield a

@pytest.fixture
def trait_lexsort_indices(trait_object):
    out = numpy.lexsort((trait_object,))
    yield out

@pytest.fixture
def mat(mat_float64, taxa_object, taxa_grp_int64, trait_object):
    out = DenseSquareTaxaSquareTraitMatrix(
        mat = mat_float64,
        taxa = taxa_object,
        taxa_grp = taxa_grp_int64,
        trait = trait_object,
    )
    out.group_taxa()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseSquareTaxaSquareTraitMatrix)

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_nsquare_is_concrete():
    assert_concrete_property(DenseSquareTaxaSquareTraitMatrix, "nsquare")

def test_square_axes_is_concrete():
    assert_concrete_property(DenseSquareTaxaSquareTraitMatrix, "square_axes")

def test_square_axes_len_is_concrete():
    assert_concrete_property(DenseSquareTaxaSquareTraitMatrix, "square_axes_len")

def test_square_taxa_axes_is_concrete():
    assert_concrete_property(DenseSquareTaxaSquareTraitMatrix, "square_taxa_axes")

def test_square_trait_axes_is_concrete():
    assert_concrete_property(DenseSquareTaxaSquareTraitMatrix, "square_trait_axes")

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "__init__")

def test_is_square_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "is_square")

def test_adjoin_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "adjoin")

def test_delete_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "delete")

def test_insert_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "insert")

def test_select_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "select")

def test_concat_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "concat")

def test_append_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "append")

def test_remove_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "remove")

def test_incorp_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "incorp")

def test_lexsort_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "lexsort")

def test_reorder_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "reorder")

def test_sort_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "sort")

def test_group_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "group")

def test_ungroup_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "ungroup")

def test_is_grouped_is_concrete():
    assert_concrete_method(DenseSquareTaxaSquareTraitMatrix, "is_grouped")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########### Matrix element copy-on-manipulation ############
def test_adjoin_cls(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test adjoin taxa
    m = mat.adjoin(mat, axis = mat.taxa_axis)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, mat.taxa_axis)
    tmp2 = numpy.append(b, a, mat.taxa_axis)
    mtrue = numpy.append(tmp1, tmp2, mat.taxa_axis+1)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    # test adjoin trait
    m = mat.adjoin(mat, axis = mat.trait_axis)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, mat.trait_axis)
    tmp2 = numpy.append(b, a, mat.trait_axis)
    mtrue = numpy.append(tmp1, tmp2, mat.trait_axis+1)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_adjoin_taxa_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test adjoin taxa
    m = mat.adjoin(mat_float64, axis = mat.taxa_axis, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, mat.taxa_axis)
    tmp2 = numpy.append(b, a, mat.taxa_axis)
    mtrue = numpy.append(tmp1, tmp2, mat.taxa_axis+1)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    # test adjoin trait
    m = mat.adjoin(mat_float64, axis = mat.trait_axis, trait = trait_object)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, mat.trait_axis)
    tmp2 = numpy.append(b, a, mat.trait_axis)
    mtrue = numpy.append(tmp1, tmp2, mat.trait_axis+1)
    assert numpy.array_equal(m.mat, mtrue, equal_nan = True)
    assert numpy.all(m.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_delete_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test delete taxa
    obj = slice(0,2,None)
    m = mat.delete(obj, axis = mat.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test delete trait
    obj = slice(0,2,None)
    m = mat.delete(obj, axis = mat.trait_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test delete taxa
    obj = 0
    m = mat.delete(obj, axis = mat.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test delete trait
    obj = 0
    m = mat.delete(obj, axis = mat.trait_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

def test_delete_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test delete taxa
    obj = [0,1,2]
    m = mat.delete(obj, axis = mat.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test delete trait
    obj = [0,1,2]
    m = mat.delete(obj, axis = mat.trait_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.delete(trait_object, obj, axis = 0))

# TODO: these tests
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

def test_select_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test select taxa
    obj = [0,0,1]
    m = mat.select(obj, axis = mat.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_taxa_axes:
        mtrue = numpy.take(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.taxa == numpy.take(taxa_object, obj, axis = 0))
    assert numpy.all(m.taxa_grp == numpy.take(taxa_grp_int64, obj, axis = 0))
    # test select trait
    obj = [0,0,1]
    m = mat.select(obj, axis = mat.trait_axis)
    mtrue = mat_float64.copy()
    for axis in mat.square_trait_axes:
        mtrue = numpy.take(mtrue, obj, axis = axis)
    assert numpy.all(m.mat == mtrue)
    assert numpy.all(m.trait == numpy.take(trait_object, obj, axis = 0))

# def test_concat_taxa_cls(mat, mat_float64, taxa_object, taxa_grp_int64):
#     obj = [mat, mat]
#     m = mat.concat_taxa(obj)
#     assert numpy.all(m.mat == numpy.concatenate([mat_float64,mat_float64], axis = mat.taxa_axis))
#     assert numpy.all(m.taxa == numpy.concatenate([taxa_object, taxa_object], axis = 0))
#     assert numpy.all(m.taxa_grp == numpy.concatenate([taxa_grp_int64, taxa_grp_int64], axis = 0))

########### Matrix element in-place-manipulation ###########
def test_append_cls(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test append taxa
    tmp = mat.deepcopy()
    tmp.append(tmp, axis = tmp.taxa_axis)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, tmp.taxa_axis)
    tmp2 = numpy.append(b, a, tmp.taxa_axis)
    mtrue = numpy.append(tmp1, tmp2, tmp.taxa_axis+1)
    assert numpy.array_equal(tmp.mat, mtrue, equal_nan = True)
    assert numpy.all(tmp.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(tmp.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    # test append trait
    tmp = mat.deepcopy()
    tmp.append(tmp, axis = tmp.trait_axis)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, tmp.trait_axis)
    tmp2 = numpy.append(b, a, tmp.trait_axis)
    mtrue = numpy.append(tmp1, tmp2, tmp.trait_axis+1)
    assert numpy.array_equal(tmp.mat, mtrue, equal_nan = True)
    print(mat.trait)
    print(tmp.trait)
    assert numpy.all(tmp.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_append_ndarray(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test append taxa
    tmp = mat.deepcopy()
    tmp.append(mat_float64, axis = tmp.taxa_axis, taxa = taxa_object, taxa_grp = taxa_grp_int64)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, tmp.taxa_axis)
    tmp2 = numpy.append(b, a, tmp.taxa_axis)
    mtrue = numpy.append(tmp1, tmp2, tmp.taxa_axis+1)
    assert numpy.array_equal(tmp.mat, mtrue, equal_nan = True)
    assert numpy.all(tmp.taxa == numpy.append(taxa_object, taxa_object, axis = 0))
    assert numpy.all(tmp.taxa_grp == numpy.append(taxa_grp_int64, taxa_grp_int64, axis = 0))
    # test append trait
    tmp = mat.deepcopy()
    tmp.append(mat_float64, axis = tmp.trait_axis, trait = trait_object)
    a = mat_float64.copy()
    b = mat_float64.copy()
    b.fill(numpy.nan)
    tmp1 = numpy.append(a, b, tmp.trait_axis)
    tmp2 = numpy.append(b, a, tmp.trait_axis)
    mtrue = numpy.append(tmp1, tmp2, tmp.trait_axis+1)
    assert numpy.array_equal(tmp.mat, mtrue, equal_nan = True)
    assert numpy.all(tmp.trait == numpy.append(trait_object, trait_object, axis = 0))

def test_remove_cls_slice(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test remove taxa
    obj = slice(0,2,None)
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(tmp.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test remove trait
    obj = slice(0,2,None)
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.trait_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_cls_int(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test remove taxa
    obj = 0
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(tmp.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test remove trait
    obj = 0
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.trait_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.trait == numpy.delete(trait_object, obj, axis = 0))

def test_remove_cls_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, trait_object):
    # test remove taxa
    obj = [0,1,2]
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.taxa_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_taxa_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.taxa == numpy.delete(taxa_object, obj, axis = 0))
    assert numpy.all(tmp.taxa_grp == numpy.delete(taxa_grp_int64, obj, axis = 0))
    # test remove trait
    obj = [0,1,2]
    tmp = mat.deepcopy()
    tmp.remove(obj, axis = tmp.trait_axis)
    mtrue = mat_float64.copy()
    for axis in tmp.square_trait_axes:
        mtrue = numpy.delete(mtrue, obj, axis = axis)
    assert numpy.all(tmp.mat == mtrue)
    assert numpy.all(tmp.trait == numpy.delete(trait_object, obj, axis = 0))

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
def test_lexsort_None(mat, taxa_lexsort_indices, trait_lexsort_indices):
    # test lexsort taxa
    ix = mat.lexsort(keys = None, axis = mat.taxa_axis)
    assert numpy.all(ix == taxa_lexsort_indices)
    # test lexsort trait
    ix = mat.lexsort(keys = None, axis = mat.trait_axis)
    assert numpy.all(ix == trait_lexsort_indices)

def test_lexsort_tuple(mat, taxa_object, taxa_grp_int64, taxa_lexsort_indices, trait_object, trait_lexsort_indices):
    # test lexsort taxa
    ix = mat.lexsort_taxa(keys = (taxa_object, taxa_grp_int64))
    assert numpy.all(ix == taxa_lexsort_indices)
    # test lexsort trait
    ix = mat.lexsort_trait(keys = (trait_object,))
    assert numpy.all(ix == trait_lexsort_indices)

def test_reorder_array_like(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices, trait_object, trait_lexsort_indices):
    # test reorder taxa
    mat.reorder(taxa_lexsort_indices, axis = mat.taxa_axis)
    tmp = mat_float64[taxa_lexsort_indices,:,:,:]
    tmp = tmp[:,taxa_lexsort_indices,:,:]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])
    # test reorder trait
    mat.reorder(trait_lexsort_indices, axis = mat.trait_axis)
    tmp = mat_float64[:,:,trait_lexsort_indices,:]
    tmp = tmp[:,:,:,trait_lexsort_indices]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_sort_None(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices, trait_object, trait_lexsort_indices):
    # test sort taxa
    mat.sort(keys = None, axis = mat.taxa_axis)
    tmp = mat_float64[taxa_lexsort_indices,:,:,:]
    tmp = tmp[:,taxa_lexsort_indices,:,:]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])
    # test sort trait
    mat.sort(keys = None, axis = mat.taxa_axis)
    tmp = mat_float64[:,:,trait_lexsort_indices,:]
    tmp = tmp[:,:,:,trait_lexsort_indices]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_sort_tuple(mat, mat_float64, taxa_object, taxa_grp_int64, taxa_lexsort_indices, trait_object, trait_lexsort_indices):
    # test sort taxa
    mat.sort(keys = (taxa_object, taxa_grp_int64), axis = mat.taxa_axis)
    tmp = mat_float64[taxa_lexsort_indices,:,:,:]
    tmp = tmp[:,taxa_lexsort_indices,:,:]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.taxa == taxa_object[taxa_lexsort_indices])
    assert numpy.all(mat.taxa_grp == taxa_grp_int64[taxa_lexsort_indices])
    # test sort trait
    mat.sort(keys = (trait_object,), axis = mat.trait_axis)
    tmp = mat_float64[:,:,trait_lexsort_indices,:]
    tmp = tmp[:,:,:,trait_lexsort_indices]
    assert numpy.all(mat.mat == tmp)
    assert numpy.all(mat.trait == trait_object[trait_lexsort_indices])

def test_group(mat, taxa_grp_name_int64, taxa_grp_stix_int64, taxa_grp_spix_int64, taxa_grp_len_int64):
    # test group taxa
    mat.group(axis = mat.taxa_axis)
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)
    # test group trait
    with pytest.raises(ValueError):
        mat.group(axis = mat.trait_axis)

def test_ungroup(mat, taxa_grp_name_int64, taxa_grp_stix_int64, taxa_grp_spix_int64, taxa_grp_len_int64):
    # test ungroup taxa
    mat.group(axis = mat.taxa_axis)
    assert numpy.all(mat.taxa_grp_name == taxa_grp_name_int64)
    assert numpy.all(mat.taxa_grp_stix == taxa_grp_stix_int64)
    assert numpy.all(mat.taxa_grp_spix == taxa_grp_spix_int64)
    assert numpy.all(mat.taxa_grp_len == taxa_grp_len_int64)
    mat.ungroup(axis = mat.taxa_axis)
    assert mat.taxa_grp_name is None
    assert mat.taxa_grp_stix is None
    assert mat.taxa_grp_spix is None
    assert mat.taxa_grp_len is None
    # test ungroup trait
    with pytest.raises(ValueError):
        mat.ungroup(axis = mat.trait_axis)

def test_is_grouped(mat):
    # test is_grouped taxa
    mat.group(axis = mat.taxa_axis)
    assert mat.is_grouped_taxa() == (
        (mat.taxa_grp_name is not None) and
        (mat.taxa_grp_stix is not None) and
        (mat.taxa_grp_spix is not None) and
        (mat.taxa_grp_len is not None)
    )
    # test is_grouped trait
    with pytest.raises(ValueError):
        mat.is_grouped(axis = mat.trait_axis)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSquareTaxaSquareTraitMatrix_is_concrete():
    assert_concrete_function(check_is_DenseSquareTaxaSquareTraitMatrix)

def test_check_is_DenseSquareTaxaSquareTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseSquareTaxaSquareTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseSquareTaxaSquareTraitMatrix(None, "mat")
