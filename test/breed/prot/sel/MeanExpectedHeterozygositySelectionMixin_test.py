import numpy
import pytest
from pybrops.breed.prot.sel.MeanExpectedHeterozygositySelection import MeanExpectedHeterozygositySelectionMixin
from pybrops.test.assert_python import assert_property_isconcrete, assert_class_documentation, assert_class_ismixin, not_raises
from .common_fixtures_large import *

################################ Test fixtures #################################

@pytest.fixture
def selmix():
    out = MeanExpectedHeterozygositySelectionMixin()
    yield out

################### Test class abstract/concrete properties ####################
def test_MeanExpectedHeterozygositySelectionMixin_is_mixin():
    assert_class_ismixin(MeanExpectedHeterozygositySelectionMixin)

############################## Test class docstring ############################
def test_MeanExpectedHeterozygositySelectionMixin_docstring():
    assert_class_documentation(MeanExpectedHeterozygositySelectionMixin)

############################ Test class properties #############################

############################## Test class methods ##############################
def test_init(selmix):
    assert selmix is not None
