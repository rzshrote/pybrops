import numpy
import pytest
from pybrops.breed.prot.sel.OptimalHaploidValueSelection import OptimalHaploidValueSelectionMixin
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_nhaploblk,
        common_unique_parents
    ):
    out = OptimalHaploidValueSelectionMixin()
    out.ntrait = common_ntrait
    out.nhaploblk = common_nhaploblk
    out.unique_parents = common_unique_parents
    yield out

################### Test class abstract/concrete properties ####################
def test_OptimalHaploidValueSelectionMixin_is_mixin():
    assert_class_ismixin(OptimalHaploidValueSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(OptimalHaploidValueSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_OptimalHaploidValueSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(OptimalHaploidValueSelectionMixin, "ntrait")

def test_ntrait_fget(selmix, common_ntrait):
    assert selmix.ntrait == common_ntrait

def test_ntrait_fset(selmix, common_ntrait):
    with not_raises(Exception):
        selmix.ntrait = int(common_ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int8(common_ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int16(common_ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int32(common_ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int64(common_ntrait)

def test_ntrait_fset_TypeError(selmix, common_ntrait):
    with pytest.raises(TypeError):
        selmix.ntrait = object()
    with pytest.raises(TypeError):
        selmix.ntrait = None
    with pytest.raises(TypeError):
        selmix.ntrait = str(common_ntrait)
    with pytest.raises(TypeError):
        selmix.ntrait = numpy.array([common_ntrait], dtype=int)

def test_ntrait_fset_ValueError(selmix, common_ntrait):
    with pytest.raises(ValueError):
        selmix.ntrait = -common_ntrait
    with pytest.raises(ValueError):
        selmix.ntrait = int(0)

def test_ntrait_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.ntrait

### nhaploblk ###
def test_OptimalHaploidValueSelectionMixin_nhaploblk_is_concrete():
    assert_property_isconcrete(OptimalHaploidValueSelectionMixin, "nhaploblk")

def test_nhaploblk_fget(selmix, common_nhaploblk):
    assert selmix.nhaploblk == common_nhaploblk

def test_nhaploblk_fset(selmix, common_nhaploblk):
    with not_raises(Exception):
        selmix.nhaploblk = int(common_nhaploblk)
    with not_raises(Exception):
        selmix.nhaploblk = numpy.int8(common_nhaploblk)
    with not_raises(Exception):
        selmix.nhaploblk = numpy.int16(common_nhaploblk)
    with not_raises(Exception):
        selmix.nhaploblk = numpy.int32(common_nhaploblk)
    with not_raises(Exception):
        selmix.nhaploblk = numpy.int64(common_nhaploblk)

def test_nhaploblk_fset_TypeError(selmix, common_nhaploblk):
    with pytest.raises(TypeError):
        selmix.nhaploblk = object()
    with pytest.raises(TypeError):
        selmix.nhaploblk = None
    with pytest.raises(TypeError):
        selmix.nhaploblk = str(common_nhaploblk)
    with pytest.raises(TypeError):
        selmix.nhaploblk = numpy.array([common_nhaploblk], dtype=int)

def test_nhaploblk_fset_ValueError(selmix, common_nhaploblk):
    with pytest.raises(ValueError):
        selmix.nhaploblk = -common_nhaploblk
    with pytest.raises(ValueError):
        selmix.nhaploblk = int(0)

def test_nhaploblk_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.nhaploblk

### unique_parents ###
def test_OptimalHaploidValueSelectionMixin_unique_parents_is_concrete():
    assert_property_isconcrete(OptimalHaploidValueSelectionMixin, "unique_parents")

def test_unique_parents_fget(selmix, common_unique_parents):
    assert isinstance(selmix.unique_parents, bool)
    assert selmix.unique_parents == common_unique_parents

def test_unique_parents_fset(selmix, common_unique_parents):
    with not_raises(Exception):
        selmix.unique_parents = common_unique_parents
    with not_raises(Exception):
        selmix.unique_parents = True
    with not_raises(Exception):
        selmix.unique_parents = False

def test_unique_parents_fset_TypeError(selmix, common_unique_parents):
    with pytest.raises(TypeError):
        selmix.unique_parents = object()
    with pytest.raises(TypeError):
        selmix.unique_parents = None
    with pytest.raises(TypeError):
        selmix.unique_parents = str(common_unique_parents)
    with pytest.raises(TypeError):
        selmix.unique_parents = numpy.array([common_unique_parents], dtype=bool)

def test_unique_parents_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.unique_parents

