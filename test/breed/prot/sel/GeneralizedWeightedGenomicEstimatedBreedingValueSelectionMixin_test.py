import numpy
import pytest
from pybrops.breed.prot.sel.GeneralizedWeightedGenomicEstimatedBreedingValueSelection import GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_alpha
    ):
    out = GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin()
    out.ntrait = common_ntrait
    out.alpha = common_alpha
    yield out

################### Test class abstract/concrete properties ####################
def test_GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin_is_mixin():
    assert_mixin_class(GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin_ntrait_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin, "ntrait")

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

### alpha ###
def test_GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin_alpha_is_concrete():
    assert_concrete_property(GeneralizedWeightedGenomicEstimatedBreedingValueSelectionMixin, "alpha")

def test_alpha_fget(selmix, common_alpha):
    assert selmix.alpha == common_alpha

def test_alpha_fset(selmix, common_alpha):
    with not_raises(Exception):
        selmix.alpha = float(common_alpha)
    with not_raises(Exception):
        selmix.alpha = numpy.float32(common_alpha)
    with not_raises(Exception):
        selmix.alpha = numpy.float64(common_alpha)

def test_alpha_fset_TypeError(selmix, common_alpha):
    with pytest.raises(TypeError):
        selmix.alpha = object()
    with pytest.raises(TypeError):
        selmix.alpha = None
    with pytest.raises(TypeError):
        selmix.alpha = str(common_alpha)
    with pytest.raises(TypeError):
        selmix.alpha = numpy.array([common_alpha], dtype=int)

def test_alpha_fset_ValueError(selmix, common_alpha):
    with pytest.raises(ValueError):
        selmix.alpha = -common_alpha
    with pytest.raises(ValueError):
        selmix.alpha = float(-2.718)
    with pytest.raises(ValueError):
        selmix.alpha = float(2.718)

def test_alpha_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.alpha
