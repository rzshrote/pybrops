import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from .common_fixtures import *

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
    assert_class_documentation(AdditiveLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_misc_is_abstract():
    assert_property_isabstract(AdditiveLinearGenomicModel, "u_misc")

def test_u_a_is_abstract():
    assert_property_isabstract(AdditiveLinearGenomicModel, "u_a")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_AdditiveLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_AdditiveLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_AdditiveLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_AdditiveLinearGenomicModel(None, "gmod")
