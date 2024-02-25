import numpy
import pytest
from pybrops.breed.prot.sel.ExpectedMaximumBreedingValueSelection import ExpectedMaximumBreedingValueSelectionMixin
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_nrep,
        common_mateprot,
        common_unique_parents
    ):
    out = ExpectedMaximumBreedingValueSelectionMixin()
    out.ntrait = common_ntrait
    out.nrep = common_nrep
    out.mateprot = common_mateprot
    out.unique_parents = common_unique_parents
    yield out

################### Test class abstract/concrete properties ####################
def test_ExpectedMaximumBreedingValueSelectionMixin_is_mixin():
    assert_class_ismixin(ExpectedMaximumBreedingValueSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(ExpectedMaximumBreedingValueSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_ExpectedMaximumBreedingValueSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueSelectionMixin, "ntrait")

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

### nrep ###
def test_ExpectedMaximumBreedingValueSelectionMixin_nrep_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueSelectionMixin, "nrep")

def test_nrep_fget(selmix, common_nrep):
    assert selmix.nrep == common_nrep

def test_nrep_fset(selmix, common_nrep):
    with not_raises(Exception):
        selmix.nrep = int(common_nrep)
    with not_raises(Exception):
        selmix.nrep = numpy.int8(common_nrep)
    with not_raises(Exception):
        selmix.nrep = numpy.int16(common_nrep)
    with not_raises(Exception):
        selmix.nrep = numpy.int32(common_nrep)
    with not_raises(Exception):
        selmix.nrep = numpy.int64(common_nrep)

def test_nrep_fset_TypeError(selmix, common_nrep):
    with pytest.raises(TypeError):
        selmix.nrep = object()
    with pytest.raises(TypeError):
        selmix.nrep = None
    with pytest.raises(TypeError):
        selmix.nrep = str(common_nrep)
    with pytest.raises(TypeError):
        selmix.nrep = numpy.array([common_nrep], dtype=int)

def test_nrep_fset_ValueError(selmix, common_nrep):
    with pytest.raises(ValueError):
        selmix.nrep = -common_nrep
    with pytest.raises(ValueError):
        selmix.nrep = int(0)

def test_nrep_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.nrep

### mateprot ###
def test_ExpectedMaximumBreedingValueSelectionMixin_mateprot_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueSelectionMixin, "mateprot")

def test_mateprot_fget(selmix, common_mateprot):
    assert selmix.mateprot == common_mateprot

def test_mateprot_fset(selmix, common_mateprot):
    with not_raises(Exception):
        selmix.mateprot = common_mateprot

def test_mateprot_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.mateprot = object()
    with pytest.raises(TypeError):
        selmix.mateprot = None

def test_mateprot_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.mateprot

### unique_parents ###
def test_ExpectedMaximumBreedingValueSelectionMixin_unique_parents_is_concrete():
    assert_property_isconcrete(ExpectedMaximumBreedingValueSelectionMixin, "unique_parents")

def test_unique_parents_fget(selmix, common_unique_parents):
    assert selmix.unique_parents == common_unique_parents

def test_unique_parents_fset(selmix):
    with not_raises(Exception):
        selmix.unique_parents = True
    with not_raises(Exception):
        selmix.unique_parents = False

def test_unique_parents_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.unique_parents = object()
    with pytest.raises(TypeError):
        selmix.unique_parents = None

def test_unique_parents_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.unique_parents
