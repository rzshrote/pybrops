import numpy
import pytest
from pybrops.breed.prot.sel.L2NormGenomicSelection import L2NormGenomicSelectionMixin
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix(
        common_cmatfcty
    ):
    out = L2NormGenomicSelectionMixin()
    out.cmatfcty = common_cmatfcty
    yield out

################### Test class abstract/concrete properties ####################
def test_L2NormGenomicSelectionMixin_is_mixin():
    assert_mixin_class(L2NormGenomicSelectionMixin)

############################## Test class docstring ############################
def test_L2NormGenomicSelectionMixin_docstring():
    assert_docstring(L2NormGenomicSelectionMixin)

############################ Test class properties #############################

### cmatfcty ###
def test_L2NormGenomicSelectionMixin_cmatfcty_is_concrete():
    assert_concrete_property(L2NormGenomicSelectionMixin, "cmatfcty")

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
