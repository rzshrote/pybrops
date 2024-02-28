import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import check_is_AdditiveDominanceLinearGenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyAdditiveDominanceLinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveDominanceLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_d_is_abstract():
    assert_property_isabstract(AdditiveDominanceLinearGenomicModel, "u_d")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_AdditiveDominanceLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_AdditiveDominanceLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveDominanceLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_AdditiveDominanceLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_AdditiveDominanceLinearGenomicModel(None, "gmod")
