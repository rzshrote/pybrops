import numpy
import pytest
from pybrops.breed.prot.sel.OptimalPopulationValueSelection import OptimalPopulationValueSelectionMixin
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_nhaploblk
    ):
    out = OptimalPopulationValueSelectionMixin()
    out.ntrait = common_ntrait
    out.nhaploblk = common_nhaploblk
    yield out

################### Test class abstract/concrete properties ####################
def test_OptimalPopulationValueSelectionMixin_is_mixin():
    assert_class_ismixin(OptimalPopulationValueSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(OptimalPopulationValueSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_OptimalPopulationValueSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(OptimalPopulationValueSelectionMixin, "ntrait")

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
def test_OptimalPopulationValueSelectionMixin_nhaploblk_is_concrete():
    assert_property_isconcrete(OptimalPopulationValueSelectionMixin, "nhaploblk")

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
