import numpy
import pytest
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSelectionMixin
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_cmatfcty,
        common_unscale
    ):
    out = OptimalContributionSelectionMixin()
    out.ntrait = common_ntrait
    out.cmatfcty = common_cmatfcty
    out.unscale = common_unscale
    yield out

################### Test class abstract/concrete properties ####################
def test_OptimalContributionSelectionMixin_is_mixin():
    assert_class_ismixin(OptimalContributionSelectionMixin)

############################## Test class docstring ############################
def test_OptimalContributionSelectionMixin_docstring():
    assert_class_documentation(OptimalContributionSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_OptimalContributionSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(OptimalContributionSelectionMixin, "ntrait")

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

### cmatfcty ###
def test_OptimalContributionSelectionMixin_cmatfcty_is_concrete():
    assert_property_isconcrete(OptimalContributionSelectionMixin, "cmatfcty")

def test_cmatfcty_fget(selmix, common_cmatfcty):
    assert isinstance(selmix.cmatfcty, CoancestryMatrixFactory)
    assert selmix.cmatfcty == common_cmatfcty

def test_cmatfcty_fset(selmix, common_cmatfcty):
    with not_raises(Exception):
        selmix.cmatfcty = common_cmatfcty

def test_cmatfcty_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.cmatfcty = object()
    with pytest.raises(TypeError):
        selmix.cmatfcty = None
    with pytest.raises(TypeError):
        selmix.cmatfcty = "string"
    with pytest.raises(TypeError):
        selmix.cmatfcty = []
    with pytest.raises(TypeError):
        selmix.cmatfcty = {}
    with pytest.raises(TypeError):
        selmix.cmatfcty = numpy.array([], dtype=int)

def test_cmatfcty_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.cmatfcty

### unscale ###
def test_OptimalContributionSelectionMixin_unscale_is_concrete():
    assert_property_isconcrete(OptimalContributionSelectionMixin, "unscale")

def test_unscale_fget(selmix, common_unscale):
    assert isinstance(selmix.unscale, bool)
    assert selmix.unscale == common_unscale

def test_unscale_fset(selmix, common_unscale):
    with not_raises(Exception):
        selmix.unscale = common_unscale

def test_unscale_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.unscale = object()
    with pytest.raises(TypeError):
        selmix.unscale = None
    with pytest.raises(TypeError):
        selmix.unscale = "string"
    with pytest.raises(TypeError):
        selmix.unscale = []
    with pytest.raises(TypeError):
        selmix.unscale = {}
    with pytest.raises(TypeError):
        selmix.unscale = numpy.array([], dtype=bool)

def test_unscale_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.unscale

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
