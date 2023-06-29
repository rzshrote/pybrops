import numpy
import pytest
from pybrops.breed.prot.sel.MeanGenomicRelationshipSelection import MeanGenomicRelationshipSelectionMixin
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from pybrops.test.breed.prot.sel.common_fixtures_large import *

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

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
