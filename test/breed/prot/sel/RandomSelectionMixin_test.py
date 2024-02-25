import numpy
import pytest
from pybrops.breed.prot.sel.RandomSelection import RandomSelectionMixin
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises

################################ Test fixtures #################################

@pytest.fixture
def ntrait():
    yield 1

@pytest.fixture
def selmix(ntrait):
    out = RandomSelectionMixin()
    out.ntrait = ntrait
    yield out

################### Test class abstract/concrete properties ####################
def test_RandomSelectionMixin_is_mixin():
    assert_class_ismixin(RandomSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(RandomSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_RandomSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(RandomSelectionMixin, "ntrait")

def test_ntrait_fget(selmix, ntrait):
    assert selmix.ntrait == ntrait

def test_ntrait_fset(selmix, ntrait):
    with not_raises(Exception):
        selmix.ntrait = int(ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int8(ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int16(ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int32(ntrait)
    with not_raises(Exception):
        selmix.ntrait = numpy.int64(ntrait)

def test_ntrait_fset_TypeError(selmix, ntrait):
    with pytest.raises(TypeError):
        selmix.ntrait = object()
    with pytest.raises(TypeError):
        selmix.ntrait = None
    with pytest.raises(TypeError):
        selmix.ntrait = str(ntrait)
    with pytest.raises(TypeError):
        selmix.ntrait = numpy.array([ntrait], dtype=int)

def test_ntrait_fset_ValueError(selmix, ntrait):
    with pytest.raises(ValueError):
        selmix.ntrait = -ntrait
    with pytest.raises(ValueError):
        selmix.ntrait = int(0)

def test_ntrait_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.ntrait
