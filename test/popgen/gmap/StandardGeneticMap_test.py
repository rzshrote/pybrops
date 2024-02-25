from numbers import Integral
import numpy
import pytest
from os.path import isfile
from pybrops.popgen.gmap.StandardGeneticMap import StandardGeneticMap, check_is_StandardGeneticMap
from pybrops.test.assert_python import assert_function_isconcrete, assert_method_isconcrete, assert_property_isconcrete, assert_class_documentation, not_raises

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmap_uncorrected(shared_datadir):
    yield StandardGeneticMap.from_csv(
        shared_datadir / "McMullen_2009_US_NAM.M.egmap",
        sep='\t',
        header=0,
        vrnt_chrgrp_col="chr_grp",
        vrnt_phypos_col="chr_start",
        vrnt_genpos_col="map_pos",
        vrnt_genpos_units="M"
    )

@pytest.fixture
def gmap(shared_datadir):
    yield StandardGeneticMap.from_csv(
        shared_datadir / "McMullen_2009_US_NAM_corrected.M.egmap",
        sep='\t',
        header=0,
        vrnt_chrgrp_col="chr_grp",
        vrnt_phypos_col="chr_start",
        vrnt_genpos_col="map_pos",
        vrnt_genpos_units="M"
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(StandardGeneticMap)

################################################################################
########################### Test concrete properties ###########################
################################################################################

### nvrnt
def test_nvrnt_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "nvrnt")

def test_nvrnt_fget(gmap):
    assert isinstance(gmap.nvrnt, Integral)
    assert len(gmap) == gmap.nvrnt

def test_nvrnt_fset_AttributeError(gmap):
    with pytest.raises(AttributeError):
        gmap.nvrnt = object()

def test_nvrnt_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.nvrnt

### vrnt_chrgrp
def test_vrnt_chrgrp_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_chrgrp")

def test_vrnt_chrgrp_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_chrgrp
    assert isinstance(gmap.vrnt_chrgrp, numpy.ndarray)
    assert len(gmap) == len(gmap.vrnt_chrgrp)

def test_vrnt_chrgrp_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_chrgrp = numpy.random.randint(1, 9, len(gmap))

def test_vrnt_chrgrp_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp = object()
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp = None
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp = 1
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp = numpy.random.random(len(gmap))

def test_vrnt_chrgrp_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_chrgrp = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_chrgrp_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_chrgrp

### vrnt_phypos
def test_vrnt_phypos_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_phypos")

def test_vrnt_phypos_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_phypos
    assert isinstance(gmap.vrnt_phypos, numpy.ndarray)
    assert len(gmap) == len(gmap.vrnt_phypos)

def test_vrnt_phypos_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_phypos = numpy.random.randint(1, 9, len(gmap))

def test_vrnt_phypos_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_phypos = object()
    with pytest.raises(TypeError):
        gmap.vrnt_phypos = None
    with pytest.raises(TypeError):
        gmap.vrnt_phypos = 1
    with pytest.raises(TypeError):
        gmap.vrnt_phypos = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_phypos = numpy.random.random(len(gmap))

def test_vrnt_phypos_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_phypos = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_phypos_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_phypos

### vrnt_genpos
def test_vrnt_genpos_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_genpos")

def test_vrnt_genpos_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_genpos
    assert isinstance(gmap.vrnt_genpos, numpy.ndarray)
    assert len(gmap) == len(gmap.vrnt_genpos)

def test_vrnt_genpos_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_genpos = numpy.random.random(len(gmap))

def test_vrnt_genpos_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_genpos = object()
    with pytest.raises(TypeError):
        gmap.vrnt_genpos = None
    with pytest.raises(TypeError):
        gmap.vrnt_genpos = 1
    with pytest.raises(TypeError):
        gmap.vrnt_genpos = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_genpos = numpy.random.randint(1, 9, len(gmap))

def test_vrnt_genpos_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_genpos = numpy.random.random((len(gmap),1))

def test_vrnt_genpos_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_genpos

### vrnt_genpos_name
def test_vrnt_chrgrp_name_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_chrgrp_name")

def test_vrnt_chrgrp_name_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_chrgrp_name
    if gmap.vrnt_chrgrp_name is not None:
        assert isinstance(gmap.vrnt_chrgrp_name, numpy.ndarray)
        assert len(gmap.vrnt_chrgrp_name) == len(numpy.unique(gmap.vrnt_chrgrp))

def test_vrnt_chrgrp_name_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_chrgrp_name = numpy.random.randint(1, 9, len(gmap))
    with not_raises(Exception):
        gmap.vrnt_chrgrp_name = None

def test_vrnt_chrgrp_name_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_name = object()
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_name = 1
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_name = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_name = numpy.random.random(len(gmap))

def test_vrnt_chrgrp_name_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_chrgrp_name = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_chrgrp_name_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_chrgrp_name

### vrnt_chrgrp_stix
def test_vrnt_chrgrp_stix_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_chrgrp_stix")

def test_vrnt_chrgrp_stix_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_chrgrp_stix
    if gmap.vrnt_chrgrp_stix is not None:
        assert isinstance(gmap.vrnt_chrgrp_stix, numpy.ndarray)
        assert len(gmap.vrnt_chrgrp_stix) == len(numpy.unique(gmap.vrnt_chrgrp))

def test_vrnt_chrgrp_stix_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_chrgrp_stix = numpy.random.randint(1, 9, len(gmap))
    with not_raises(Exception):
        gmap.vrnt_chrgrp_stix = None

def test_vrnt_chrgrp_stix_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_stix = object()
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_stix = 1
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_stix = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_stix = numpy.random.random(len(gmap))

def test_vrnt_chrgrp_stix_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_chrgrp_stix = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_chrgrp_stix_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_chrgrp_stix

### vrnt_chrgrp_spix
def test_vrnt_chrgrp_spix_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_chrgrp_spix")

def test_vrnt_chrgrp_spix_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_chrgrp_spix
    if gmap.vrnt_chrgrp_spix is not None:
        assert isinstance(gmap.vrnt_chrgrp_spix, numpy.ndarray)
        assert len(gmap.vrnt_chrgrp_spix) == len(numpy.unique(gmap.vrnt_chrgrp))

def test_vrnt_chrgrp_spix_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_chrgrp_spix = numpy.random.randint(1, 9, len(gmap))
    with not_raises(Exception):
        gmap.vrnt_chrgrp_spix = None

def test_vrnt_chrgrp_spix_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_spix = object()
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_spix = 1
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_spix = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_spix = numpy.random.random(len(gmap))

def test_vrnt_chrgrp_spix_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_chrgrp_spix = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_chrgrp_spix_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_chrgrp_spix

### vrnt_chrgrp_len
def test_vrnt_chrgrp_len_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "vrnt_chrgrp_len")

def test_vrnt_chrgrp_len_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.vrnt_chrgrp_len
    if gmap.vrnt_chrgrp_len is not None:
        assert isinstance(gmap.vrnt_chrgrp_len, numpy.ndarray)
        assert len(gmap.vrnt_chrgrp_len) == len(numpy.unique(gmap.vrnt_chrgrp))

def test_vrnt_chrgrp_len_fset(gmap):
    with not_raises(Exception):
        gmap.vrnt_chrgrp_len = numpy.random.randint(1, 9, len(gmap))
    with not_raises(Exception):
        gmap.vrnt_chrgrp_len = None

def test_vrnt_chrgrp_len_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_len = object()
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_len = 1
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_len = 1.0
    with pytest.raises(TypeError):
        gmap.vrnt_chrgrp_len = numpy.random.random(len(gmap))

def test_vrnt_chrgrp_len_fset_ValueError(gmap):
    with pytest.raises(ValueError):
        gmap.vrnt_chrgrp_len = numpy.random.randint(1, 9, (len(gmap),1))

def test_vrnt_chrgrp_len_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.vrnt_chrgrp_len

### spline
def test_spline_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "spline")

def test_spline_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.spline
    if gmap.spline is not None:
        assert isinstance(gmap.spline, dict)
        assert len(gmap.spline) == len(numpy.unique(gmap.vrnt_chrgrp))

def test_spline_fset(gmap):
    with not_raises(Exception):
        gmap.spline = None
    with not_raises(Exception):
        gmap.spline = {}

def test_spline_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.spline = object()
    with pytest.raises(TypeError):
        gmap.spline = []
    with pytest.raises(TypeError):
        gmap.spline = ()
    with pytest.raises(TypeError):
        gmap.spline = "string"
    with pytest.raises(TypeError):
        gmap.spline = 1

def test_spline_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.spline

### spline_kind
def test_spline_kind_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "spline_kind")

def test_spline_kind_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.spline_kind
    if gmap.spline_kind is not None:
        assert isinstance(gmap.spline_kind, str)

def test_spline_kind_fset(gmap):
    with not_raises(Exception):
        gmap.spline_kind = None
    with not_raises(Exception):
        gmap.spline_kind = "cubic"

def test_spline_kind_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.spline_kind = object()

def test_spline_kind_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.spline_kind

### spline_fill_value
def test_spline_fill_value_is_concrete():
    assert_property_isconcrete(StandardGeneticMap, "spline_fill_value")

def test_spline_fill_value_fget(gmap):
    with not_raises(Exception):
        tmp = gmap.spline_fill_value
    if gmap.spline_fill_value is not None:
        assert isinstance(gmap.spline_fill_value, (str,numpy.ndarray))

def test_spline_fill_value_fset(gmap):
    with not_raises(Exception):
        gmap.spline_fill_value = None
    with not_raises(Exception):
        gmap.spline_fill_value = "extrapolate"
    with not_raises(Exception):
        gmap.spline_fill_value = numpy.random.random(len(gmap))

def test_spline_fill_value_fset_TypeError(gmap):
    with pytest.raises(TypeError):
        gmap.spline_fill_value = object()

def test_spline_fill_value_fdel_AttributeError(gmap):
    with pytest.raises(AttributeError):
        del gmap.spline_fill_value

################################################################################
############################# Test concrete methods ############################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "__init__")

### __len__
def test___len___is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "__len__")

def test___len__(gmap):
    assert len(gmap) == len(gmap.vrnt_chrgrp)

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "__copy__")

def test___copy__(gmap):
    tmp = gmap.__copy__()
    assert numpy.all(gmap.vrnt_chrgrp == tmp.vrnt_chrgrp)
    assert numpy.all(gmap.vrnt_phypos == tmp.vrnt_phypos)
    assert numpy.all(gmap.vrnt_genpos == tmp.vrnt_genpos)
    assert numpy.all(gmap.vrnt_chrgrp_name == tmp.vrnt_chrgrp_name)
    assert numpy.all(gmap.vrnt_chrgrp_stix == tmp.vrnt_chrgrp_stix)
    assert numpy.all(gmap.vrnt_chrgrp_spix == tmp.vrnt_chrgrp_spix)
    assert numpy.all(gmap.vrnt_chrgrp_len == tmp.vrnt_chrgrp_len)
    assert gmap.spline == tmp.spline
    assert gmap.spline_kind == tmp.spline_kind
    assert numpy.all(gmap.spline_fill_value == tmp.spline_fill_value)

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "__deepcopy__")

def test___deepcopy__(gmap):
    tmp = gmap.__deepcopy__()
    assert numpy.all(gmap.vrnt_chrgrp == tmp.vrnt_chrgrp)
    assert numpy.all(gmap.vrnt_phypos == tmp.vrnt_phypos)
    assert numpy.all(gmap.vrnt_genpos == tmp.vrnt_genpos)
    assert numpy.all(gmap.vrnt_chrgrp_name == tmp.vrnt_chrgrp_name)
    assert numpy.all(gmap.vrnt_chrgrp_stix == tmp.vrnt_chrgrp_stix)
    assert numpy.all(gmap.vrnt_chrgrp_spix == tmp.vrnt_chrgrp_spix)
    assert numpy.all(gmap.vrnt_chrgrp_len == tmp.vrnt_chrgrp_len)
    assert gmap.spline != tmp.spline # all spline objects should be different
    assert gmap.spline_kind == tmp.spline_kind
    assert numpy.all(gmap.spline_fill_value == tmp.spline_fill_value)

### copy()
def test_copy_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "copy")

def test_copy(gmap):
    tmp = gmap.copy()
    assert numpy.all(gmap.vrnt_chrgrp == tmp.vrnt_chrgrp)
    assert numpy.all(gmap.vrnt_phypos == tmp.vrnt_phypos)
    assert numpy.all(gmap.vrnt_genpos == tmp.vrnt_genpos)
    assert numpy.all(gmap.vrnt_chrgrp_name == tmp.vrnt_chrgrp_name)
    assert numpy.all(gmap.vrnt_chrgrp_stix == tmp.vrnt_chrgrp_stix)
    assert numpy.all(gmap.vrnt_chrgrp_spix == tmp.vrnt_chrgrp_spix)
    assert numpy.all(gmap.vrnt_chrgrp_len == tmp.vrnt_chrgrp_len)
    assert gmap.spline == tmp.spline
    assert gmap.spline_kind == tmp.spline_kind
    assert numpy.all(gmap.spline_fill_value == tmp.spline_fill_value)

### deepcopy()
def test_deepcopy_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "deepcopy")

def test_deepcopy(gmap):
    tmp = gmap.deepcopy()
    assert numpy.all(gmap.vrnt_chrgrp == tmp.vrnt_chrgrp)
    assert numpy.all(gmap.vrnt_phypos == tmp.vrnt_phypos)
    assert numpy.all(gmap.vrnt_genpos == tmp.vrnt_genpos)
    assert numpy.all(gmap.vrnt_chrgrp_name == tmp.vrnt_chrgrp_name)
    assert numpy.all(gmap.vrnt_chrgrp_stix == tmp.vrnt_chrgrp_stix)
    assert numpy.all(gmap.vrnt_chrgrp_spix == tmp.vrnt_chrgrp_spix)
    assert numpy.all(gmap.vrnt_chrgrp_len == tmp.vrnt_chrgrp_len)
    assert gmap.spline != tmp.spline # all spline objects should be different
    assert gmap.spline_kind == tmp.spline_kind
    assert numpy.all(gmap.spline_fill_value == tmp.spline_fill_value)

### lexsort()
def test_lexsort_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "lexsort")

def test_lexsort_None(gmap):
    output = gmap.lexsort(keys = None)
    expected = numpy.lexsort(keys = (gmap.vrnt_genpos, gmap.vrnt_phypos, gmap.vrnt_chrgrp))
    assert numpy.all(output == expected)

def test_lexsort_tuple(gmap):
    key1 = numpy.random.randint(1, 9, len(gmap))
    key2 = numpy.random.random(len(gmap))
    output = gmap.lexsort(keys = (key2, key1))
    expected = numpy.lexsort(keys = (key2, key1))
    assert numpy.all(output == expected)

### reorder()
def test_reorder_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "reorder")

def test_reorder(gmap):
    tmp = gmap.deepcopy()
    indices = numpy.arange(len(gmap))
    numpy.random.shuffle(indices)
    tmp.reorder(indices)
    assert numpy.all(tmp.vrnt_chrgrp == gmap.vrnt_chrgrp[indices])
    assert numpy.all(tmp.vrnt_phypos == gmap.vrnt_phypos[indices])
    assert numpy.all(tmp.vrnt_genpos == gmap.vrnt_genpos[indices])
    assert tmp.vrnt_chrgrp_name is None
    assert tmp.vrnt_chrgrp_stix is None
    assert tmp.vrnt_chrgrp_spix is None
    assert tmp.vrnt_chrgrp_len is None

### sort()
def test_sort_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "sort")

def test_sort(gmap):
    tmp = gmap.deepcopy()
    key = numpy.random.random(len(gmap))
    indices = numpy.lexsort((key,))
    tmp.sort(key)
    assert numpy.all(tmp.vrnt_chrgrp == gmap.vrnt_chrgrp[indices])
    assert numpy.all(tmp.vrnt_phypos == gmap.vrnt_phypos[indices])
    assert numpy.all(tmp.vrnt_genpos == gmap.vrnt_genpos[indices])
    assert tmp.vrnt_chrgrp_name is None
    assert tmp.vrnt_chrgrp_stix is None
    assert tmp.vrnt_chrgrp_spix is None
    assert tmp.vrnt_chrgrp_len is None

### group()
def test_group_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "group")

def test_group(gmap):
    tmp = gmap.deepcopy()
    indices = numpy.arange(len(gmap))
    numpy.random.shuffle(indices)
    tmp.reorder(indices)
    tmp.group()
    assert numpy.all(tmp.vrnt_chrgrp == gmap.vrnt_chrgrp)
    assert numpy.all(tmp.vrnt_phypos == gmap.vrnt_phypos)
    assert numpy.all(tmp.vrnt_genpos == gmap.vrnt_genpos)
    assert numpy.all(tmp.vrnt_chrgrp_name == gmap.vrnt_chrgrp_name)
    assert numpy.all(tmp.vrnt_chrgrp_stix == gmap.vrnt_chrgrp_stix)
    assert numpy.all(tmp.vrnt_chrgrp_spix == gmap.vrnt_chrgrp_spix)
    assert numpy.all(tmp.vrnt_chrgrp_len == gmap.vrnt_chrgrp_len)

### ungroup()
def test_ungroup_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "ungroup")

def test_ungroup(gmap):
    assert gmap.is_grouped()
    gmap.ungroup()
    assert not gmap.is_grouped()

### is_grouped()
def test_is_grouped_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "is_grouped")

def test_is_grouped(gmap):
    tmp = gmap.deepcopy()
    indices = numpy.arange(len(gmap))
    numpy.random.shuffle(indices)
    tmp.reorder(indices)
    assert not tmp.is_grouped()
    tmp.group()
    assert tmp.is_grouped()

### congruence()
def test_congruence_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "congruence")

def test_congruence(gmap, gmap_uncorrected):
    out = gmap.congruence()
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out)
    out = gmap_uncorrected.congruence()
    assert isinstance(out, numpy.ndarray)
    assert not numpy.all(out)

### is_congruent()
def test_is_congruent_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "is_congruent")

def test_is_congruent(gmap, gmap_uncorrected):
    assert gmap.is_congruent()
    assert not gmap_uncorrected.is_congruent()

### remove_discrepancies()
def test_remove_discrepancies_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "remove_discrepancies")

def test_remove_discrepancies(gmap_uncorrected):
    tmp = gmap_uncorrected.deepcopy()
    tmp.remove_discrepancies()
    assert len(gmap_uncorrected) != len(tmp)

### remove()
def test_remove_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "remove")

def test_remove(gmap):
    # test with grouped gmap
    tmp = gmap.deepcopy()
    indices = numpy.arange(len(gmap) // 2)
    tmp.remove(indices)
    assert numpy.all(tmp.vrnt_chrgrp == numpy.delete(gmap.vrnt_chrgrp, indices))
    assert numpy.all(tmp.vrnt_phypos == numpy.delete(gmap.vrnt_phypos, indices))
    assert numpy.all(tmp.vrnt_genpos == numpy.delete(gmap.vrnt_genpos, indices))
    assert tmp.vrnt_chrgrp_name is not None
    assert tmp.vrnt_chrgrp_stix is not None
    assert tmp.vrnt_chrgrp_spix is not None
    assert tmp.vrnt_chrgrp_len is not None
    # test with ungrouped gmap
    tmp = gmap.deepcopy()
    tmp.ungroup()
    indices = numpy.arange(len(gmap) // 2)
    tmp.remove(indices)
    assert numpy.all(tmp.vrnt_chrgrp == numpy.delete(gmap.vrnt_chrgrp, indices))
    assert numpy.all(tmp.vrnt_phypos == numpy.delete(gmap.vrnt_phypos, indices))
    assert numpy.all(tmp.vrnt_genpos == numpy.delete(gmap.vrnt_genpos, indices))
    assert tmp.vrnt_chrgrp_name is None
    assert tmp.vrnt_chrgrp_stix is None
    assert tmp.vrnt_chrgrp_spix is None
    assert tmp.vrnt_chrgrp_len is None

### select()
def test_select_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "select")

def test_select(gmap):
    # test with grouped gmap
    tmp = gmap.deepcopy()
    indices = numpy.arange(len(gmap) // 2)
    tmp.select(indices)
    assert numpy.all(tmp.vrnt_chrgrp == gmap.vrnt_chrgrp[indices])
    assert numpy.all(tmp.vrnt_phypos == gmap.vrnt_phypos[indices])
    assert numpy.all(tmp.vrnt_genpos == gmap.vrnt_genpos[indices])
    assert tmp.vrnt_chrgrp_name is not None
    assert tmp.vrnt_chrgrp_stix is not None
    assert tmp.vrnt_chrgrp_spix is not None
    assert tmp.vrnt_chrgrp_len is not None
    # test with ungrouped gmap
    tmp = gmap.deepcopy()
    tmp.ungroup()
    indices = numpy.arange(len(gmap) // 2)
    tmp.select(indices)
    assert numpy.all(tmp.vrnt_chrgrp == gmap.vrnt_chrgrp[indices])
    assert numpy.all(tmp.vrnt_phypos == gmap.vrnt_phypos[indices])
    assert numpy.all(tmp.vrnt_genpos == gmap.vrnt_genpos[indices])
    assert tmp.vrnt_chrgrp_name is None
    assert tmp.vrnt_chrgrp_stix is None
    assert tmp.vrnt_chrgrp_spix is None
    assert tmp.vrnt_chrgrp_len is None

### build_spline()
def test_build_spline_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "build_spline")

def test_build_spline(gmap):
    gmap.spline = None
    assert gmap.spline == None
    kind = "linear"
    fill_value = "extrapolate"
    gmap.build_spline(kind, fill_value)
    assert gmap.spline is not None
    assert gmap.spline_kind == kind
    assert numpy.all(gmap.spline_fill_value == fill_value)

### has_spline()
def test_has_spline_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "has_spline")

def test_has_spline(gmap):
    assert gmap.has_spline()
    gmap.spline = None
    assert not gmap.has_spline()

### interp_genpos()
def test_interp_genpos_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "interp_genpos")

def test_interp_genpos(gmap):
    genpos = gmap.interp_genpos(gmap.vrnt_chrgrp, gmap.vrnt_phypos)
    assert isinstance(genpos, numpy.ndarray)
    assert numpy.all(genpos == gmap.vrnt_genpos)

### interp_gmap()
def test_interp_gmap_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "interp_gmap")

def test_interp_gmap(gmap):
    tmp = gmap.interp_gmap(gmap.vrnt_chrgrp, gmap.vrnt_phypos)
    assert isinstance(tmp, StandardGeneticMap)
    assert tmp.has_spline() == gmap.has_spline()

### gdist1g()
def test_gdist1g_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "gdist1g")

def test_gdist1g(gmap):
    out = gmap.gdist1g(gmap.vrnt_chrgrp, gmap.vrnt_genpos)
    assert isinstance(out, numpy.ndarray)
    assert out.ndim == 1
    assert out[0] == numpy.inf
    assert not numpy.all(out == numpy.inf)
    assert len(out) == len(gmap)

### gdist1p()
def test_gdist1p_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "gdist1p")

def test_gdist1p(gmap):
    out = gmap.gdist1p(gmap.vrnt_chrgrp, gmap.vrnt_phypos)
    assert isinstance(out, numpy.ndarray)
    assert out.ndim == 1
    assert out[0] == numpy.inf
    assert not numpy.all(out == numpy.inf)
    assert len(out) == len(gmap)

### gdist2g()
def test_gdist2g_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "gdist2g")

def test_gdist2g(gmap):
    out = gmap.gdist2g(gmap.vrnt_chrgrp, gmap.vrnt_genpos)
    assert isinstance(out, numpy.ndarray)
    assert out.ndim == 2
    assert not numpy.all(out == numpy.inf)
    assert out.shape == (len(gmap),len(gmap))

### gdist2p()
def test_gdist2p_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "gdist2p")

def test_gdist2p(gmap):
    out = gmap.gdist2p(gmap.vrnt_chrgrp, gmap.vrnt_phypos)
    assert isinstance(out, numpy.ndarray)
    assert out.ndim == 2
    assert not numpy.all(out == numpy.inf)
    assert out.shape == (len(gmap),len(gmap))

### to_pandas()
def test_to_pandas_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "to_pandas")

def test_to_pandas(gmap):
    df = gmap.to_pandas("chr","pos","cM")
    assert "chr" in df.columns
    assert "pos" in df.columns
    assert "cM" in df.columns

### to_csv()
def test_to_csv_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "to_csv")

def test_to_csv(gmap):
    fout = "genetic_map.csv"
    gmap.to_csv(fout)
    assert isfile(fout)

### from_pandas()
def test_from_pandas_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "from_pandas")

def test_from_pandas(gmap):
    df = gmap.to_pandas()
    tmp = StandardGeneticMap.from_pandas(df)
    assert isinstance(tmp, StandardGeneticMap)

### from_csv()
def test_from_csv_is_concrete():
    assert_method_isconcrete(StandardGeneticMap, "from_csv")

def test_from_csv(gmap):
    fin = "genetic_map.csv"
    tmp = StandardGeneticMap.from_csv(fin)
    assert isinstance(tmp, StandardGeneticMap)

################################################################################
############################ Test concrete functions ###########################
################################################################################

### check_is_StandardGeneticMap
def test_check_is_StandardGeneticMap_is_concrete():
    assert_function_isconcrete(check_is_StandardGeneticMap)

def test_check_is_StandardGeneticMap():
    with pytest.raises(TypeError):
        check_is_StandardGeneticMap(None, "None")
