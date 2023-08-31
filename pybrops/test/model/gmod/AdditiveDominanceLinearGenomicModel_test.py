import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import check_is_AdditiveDominanceLinearGenomicModel
from pybrops.test.model.gmod.common_fixtures import *

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
    assert_docstring(AdditiveDominanceLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveDominanceLinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_u_d_is_abstract():
    assert_abstract_property(AdditiveDominanceLinearGenomicModel, "u_d")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_AdditiveDominanceLinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_AdditiveDominanceLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveDominanceLinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_AdditiveDominanceLinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_AdditiveDominanceLinearGenomicModel(None, "gmod")
