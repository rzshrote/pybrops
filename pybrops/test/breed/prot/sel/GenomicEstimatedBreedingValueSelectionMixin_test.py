import numpy
import pytest
from pybrops.breed.prot.sel.GenomicEstimatedBreedingValueSelection import GenomicEstimatedBreedingValueSelectionMixin
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from pybrops.test.breed.prot.sel.common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait
    ):
    out = GenomicEstimatedBreedingValueSelectionMixin()
    out.ntrait = common_ntrait
    yield out

################### Test class abstract/concrete properties ####################
def test_GenomicEstimatedBreedingValueSelectionMixin_is_mixin():
    assert_mixin_class(GenomicEstimatedBreedingValueSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(GenomicEstimatedBreedingValueSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_GenomicEstimatedBreedingValueSelectionMixin_ntrait_is_concrete():
    assert_concrete_property(GenomicEstimatedBreedingValueSelectionMixin, "ntrait")

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
