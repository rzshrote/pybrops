import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel
from pybrops.model.gmod.LinearGenomicModel import is_LinearGenomicModel
from pybrops.model.gmod.LinearGenomicModel import check_is_LinearGenomicModel

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield LinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(LinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(LinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_beta_is_abstract():
    generic_assert_abstract_property(LinearGenomicModel, "beta")

def test_u_is_abstract():
    generic_assert_abstract_property(LinearGenomicModel, "u")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_LinearGenomicModel_is_concrete():
    generic_assert_concrete_function(is_LinearGenomicModel)

def test_check_is_LinearGenomicModel_is_concrete():
    generic_assert_concrete_function(check_is_LinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_LinearGenomicModel(gmod):
    assert is_LinearGenomicModel(gmod)

def test_check_is_LinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_LinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_LinearGenomicModel(None, "gmod")
