import numpy
import pytest
from pybrops.breed.prot.sel.MeanGenomicRelationshipSelection import MeanGenomicRelationshipSelectionMixin
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_cmatfcty
    ):
    out = MeanGenomicRelationshipSelectionMixin()
    out.cmatfcty = common_cmatfcty
    yield out

################### Test class abstract/concrete properties ####################
def test_MeanGenomicRelationshipSelectionMixin_is_mixin():
    assert_mixin_class(MeanGenomicRelationshipSelectionMixin)

############################## Test class docstring ############################
def test_MeanGenomicRelationshipSelectionMixin_docstring():
    assert_docstring(MeanGenomicRelationshipSelectionMixin)

############################ Test class properties #############################

### cmatfcty ###
def test_MeanGenomicRelationshipSelectionMixin_cmatfcty_is_concrete():
    assert_concrete_property(MeanGenomicRelationshipSelectionMixin, "cmatfcty")

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

def test_cmatfcty_fdel(selmix):
    with pytest.raises(AttributeError):
        del selmix.cmatfcty

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
