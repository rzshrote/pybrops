import numpy
import pytest
from pybrops.breed.prot.sel.OptimalContributionSelection import OptimalContributionSelectionMixin
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from pybrops.test.breed.prot.sel.common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_cmatfcty,
        common_descale
    ):
    out = OptimalContributionSelectionMixin()
    out.ntrait = common_ntrait
    out.cmatfcty = common_cmatfcty
    out.descale = common_descale
    yield out

################### Test class abstract/concrete properties ####################
def test_OptimalContributionSelectionMixin_is_mixin():
    assert_mixin_class(OptimalContributionSelectionMixin)

############################## Test class docstring ############################
def test_OptimalContributionSelectionMixin_docstring():
    assert_docstring(OptimalContributionSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_OptimalContributionSelectionMixin_ntrait_is_concrete():
    assert_concrete_property(OptimalContributionSelectionMixin, "ntrait")

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
    assert_concrete_property(OptimalContributionSelectionMixin, "cmatfcty")

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

### descale ###
def test_OptimalContributionSelectionMixin_descale_is_concrete():
    assert_concrete_property(OptimalContributionSelectionMixin, "descale")

def test_descale_fget(selmix, common_descale):
    assert isinstance(selmix.descale, bool)
    assert selmix.descale == common_descale

def test_descale_fset(selmix, common_descale):
    with not_raises(Exception):
        selmix.descale = common_descale

def test_descale_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.descale = object()
    with pytest.raises(TypeError):
        selmix.descale = None
    with pytest.raises(TypeError):
        selmix.descale = "string"
    with pytest.raises(TypeError):
        selmix.descale = []
    with pytest.raises(TypeError):
        selmix.descale = {}
    with pytest.raises(TypeError):
        selmix.descale = numpy.array([], dtype=bool)

def test_descale_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.descale

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
