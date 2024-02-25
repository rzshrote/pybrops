from numbers import Integral, Real
import numpy
import pytest
from pybrops.breed.prot.sel.UsefulnessCriterionSelection import UsefulnessCriterionSelectionMixin
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_ntrait,
        common_nself,
        common_upper_percentile,
        common_vmatfcty,
        common_gmapfn,
        common_unique_parents
    ):
    out = UsefulnessCriterionSelectionMixin()
    out.ntrait = common_ntrait
    out.nself = common_nself
    out.upper_percentile = common_upper_percentile
    out.vmatfcty = common_vmatfcty
    out.gmapfn = common_gmapfn
    out.unique_parents = common_unique_parents
    yield out

################### Test class abstract/concrete properties ####################
def test_UsefulnessCriterionSelectionMixin_is_mixin():
    assert_class_ismixin(UsefulnessCriterionSelectionMixin)

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(UsefulnessCriterionSelectionMixin)

############################ Test class properties #############################

### ntrait ###
def test_UsefulnessCriterionSelectionMixin_ntrait_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "ntrait")

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

### nself ###
def test_UsefulnessCriterionSelectionMixin_nself_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "nself")

def test_nself_fget(selmix, common_nself):
    assert isinstance(selmix.nself, Integral)
    assert selmix.nself == common_nself

def test_nself_fset(selmix, common_nself):
    with not_raises(Exception):
        selmix.nself = common_nself
    with not_raises(Exception):
        selmix.nself = int(0)
    with not_raises(Exception):
        selmix.nself = int(1)
    with not_raises(Exception):
        selmix.nself = numpy.int8(1)
    with not_raises(Exception):
        selmix.nself = numpy.int16(1)
    with not_raises(Exception):
        selmix.nself = numpy.int32(1)
    with not_raises(Exception):
        selmix.nself = numpy.int64(1)

def test_nself_fset_TypeError(selmix, common_nself):
    with pytest.raises(TypeError):
        selmix.nself = object()
    with pytest.raises(TypeError):
        selmix.nself = None
    with pytest.raises(TypeError):
        selmix.nself = str(common_nself)
    with pytest.raises(TypeError):
        selmix.nself = numpy.array([common_nself], dtype=int)

def test_nself_fset_ValueError(selmix, common_nself):
    with pytest.raises(ValueError):
        selmix.nself = int(-1)

def test_nself_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.nself

### upper_percentile ###
def test_UsefulnessCriterionSelectionMixin_upper_percentile_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "upper_percentile")

def test_upper_percentile_fget(selmix, common_upper_percentile):
    assert isinstance(selmix.upper_percentile, Real)
    assert selmix.upper_percentile >= 0.0
    assert selmix.upper_percentile <= 1.0
    assert selmix.upper_percentile == common_upper_percentile

def test_upper_percentile_fset(selmix, common_upper_percentile):
    with not_raises(Exception):
        selmix.upper_percentile = common_upper_percentile
    with not_raises(Exception):
        selmix.upper_percentile = float(0.5)
    with not_raises(Exception):
        selmix.upper_percentile = numpy.float32(0.5)
    with not_raises(Exception):
        selmix.upper_percentile = numpy.float64(0.5)

def test_upper_percentile_fset_TypeError(selmix, common_upper_percentile):
    with pytest.raises(TypeError):
        selmix.upper_percentile = object()
    with pytest.raises(TypeError):
        selmix.upper_percentile = None
    with pytest.raises(TypeError):
        selmix.upper_percentile = str(common_upper_percentile)
    with pytest.raises(TypeError):
        selmix.upper_percentile = numpy.array([common_upper_percentile], dtype=int)

def test_upper_percentile_fset_ValueError(selmix, common_upper_percentile):
    with pytest.raises(ValueError):
        selmix.upper_percentile = int(-1)
    with pytest.raises(ValueError):
        selmix.upper_percentile = int(0)
    with pytest.raises(ValueError):
        selmix.upper_percentile = int(1)
    with pytest.raises(ValueError):
        selmix.upper_percentile = float(-1)
    with pytest.raises(ValueError):
        selmix.upper_percentile = float(0)
    with pytest.raises(ValueError):
        selmix.upper_percentile = float(1)

def test_upper_percentile_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.upper_percentile

### selection_intensity ###
def test_UsefulnessCriterionSelectionMixin_selection_intensity_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "selection_intensity")

def test_selection_intensity_fget(selmix):
    assert isinstance(selmix.selection_intensity, Real)
    assert selmix.selection_intensity >= -numpy.inf
    assert selmix.selection_intensity <= numpy.inf

def test_selection_intensity_fset_AttributeError(selmix):
    with pytest.raises(AttributeError):
        selmix.selection_intensity = float(0.5)

def test_selection_intensity_fdel_AttributeError(selmix):
    with pytest.raises(AttributeError):
        del selmix.selection_intensity

### vmatfcty ###
def test_UsefulnessCriterionSelectionMixin_vmatfcty_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "vmatfcty")

def test_vmatfcty_fget(selmix, common_vmatfcty):
    assert isinstance(selmix.vmatfcty, GeneticVarianceMatrixFactory)
    assert selmix.vmatfcty == common_vmatfcty

def test_vmatfcty_fset(selmix, common_vmatfcty):
    with not_raises(Exception):
        selmix.vmatfcty = common_vmatfcty

def test_vmatfcty_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.vmatfcty = object()
    with pytest.raises(TypeError):
        selmix.vmatfcty = None
    with pytest.raises(TypeError):
        selmix.vmatfcty = "string"

def test_vmatfcty_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.vmatfcty

### gmapfn ###
def test_UsefulnessCriterionSelectionMixin_gmapfn_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "gmapfn")

def test_gmapfn_fget(selmix, common_gmapfn):
    assert isinstance(selmix.gmapfn, GeneticMapFunction)
    assert selmix.gmapfn == common_gmapfn

def test_gmapfn_fset(selmix, common_gmapfn):
    with not_raises(Exception):
        selmix.gmapfn = common_gmapfn

def test_gmapfn_fset_TypeError(selmix):
    with pytest.raises(TypeError):
        selmix.gmapfn = object()
    with pytest.raises(TypeError):
        selmix.gmapfn = None
    with pytest.raises(TypeError):
        selmix.gmapfn = "string"

def test_gmapfn_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.gmapfn

### unique_parents ###
def test_UsefulnessCriterionSelectionMixin_unique_parents_is_concrete():
    assert_property_isconcrete(UsefulnessCriterionSelectionMixin, "unique_parents")

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

