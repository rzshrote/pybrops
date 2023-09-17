import numpy
import pytest
from pybrops.breed.prot.sel.MeanExpectedHeterozygositySelection import MeanExpectedHeterozygositySelectionMixin
from pybrops.test.assert_python import assert_concrete_property, assert_docstring, assert_mixin_class, not_raises
from pybrops.test.breed.prot.sel.common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix():
    out = MeanExpectedHeterozygositySelectionMixin()
    yield out

################### Test class abstract/concrete properties ####################
def test_MeanExpectedHeterozygositySelectionMixin_is_mixin():
    assert_mixin_class(MeanExpectedHeterozygositySelectionMixin)

############################## Test class docstring ############################
def test_MeanExpectedHeterozygositySelectionMixin_docstring():
    assert_docstring(MeanExpectedHeterozygositySelectionMixin)

############################ Test class properties #############################

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
