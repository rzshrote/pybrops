import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel
from pybrops.model.gmod.LinearGenomicModel import check_is_LinearGenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyLinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(LinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(LinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_beta_is_abstract():
    assert_abstract_property(LinearGenomicModel, "beta")

def test_u_is_abstract():
    assert_abstract_property(LinearGenomicModel, "u")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_LinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_LinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_LinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_LinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_LinearGenomicModel(None, "gmod")
