import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from pybrops.test.model.gmod.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyAdditiveLinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(AdditiveLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveLinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_misc_is_abstract():
    assert_abstract_property(AdditiveLinearGenomicModel, "u_misc")

def test_u_a_is_abstract():
    assert_abstract_property(AdditiveLinearGenomicModel, "u_a")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_AdditiveLinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_AdditiveLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_AdditiveLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_AdditiveLinearGenomicModel(None, "gmod")
