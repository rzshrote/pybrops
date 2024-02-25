import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.CoancestryLinearGenomicModel import CoancestryLinearGenomicModel
from pybrops.model.gmod.CoancestryLinearGenomicModel import check_is_CoancestryLinearGenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyCoancestryLinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(CoancestryLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(CoancestryLinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_misc_is_abstract():
    assert_property_isabstract(CoancestryLinearGenomicModel, "u_misc")

def test_u_c_is_abstract():
    assert_property_isabstract(CoancestryLinearGenomicModel, "u_c")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_CoancestryLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_CoancestryLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_CoancestryLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_CoancestryLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_CoancestryLinearGenomicModel(None, "gmod")
