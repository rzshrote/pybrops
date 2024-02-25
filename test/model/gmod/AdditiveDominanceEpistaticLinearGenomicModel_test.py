import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import AdditiveDominanceEpistaticLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import check_is_AdditiveDominanceEpistaticLinearGenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyAdditiveDominanceEpistaticLinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveDominanceEpistaticLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(AdditiveDominanceEpistaticLinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_i_is_abstract():
    assert_property_isabstract(AdditiveDominanceEpistaticLinearGenomicModel, "u_i")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_AdditiveDominanceEpistaticLinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_AdditiveDominanceEpistaticLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveDominanceEpistaticLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_AdditiveDominanceEpistaticLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_AdditiveDominanceEpistaticLinearGenomicModel(None, "gmod")
