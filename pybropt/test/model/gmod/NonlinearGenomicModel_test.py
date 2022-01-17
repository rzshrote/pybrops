import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybropt.model.gmod.NonlinearGenomicModel import is_NonlinearGenomicModel
from pybropt.model.gmod.NonlinearGenomicModel import check_is_NonlinearGenomicModel
from pybropt.model.gmod.NonlinearGenomicModel import cond_check_is_NonlinearGenomicModel

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield NonlinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(NonlinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(NonlinearGenomicModel, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_NonlinearGenomicModel_is_concrete():
    generic_assert_concrete_function(is_NonlinearGenomicModel)

def test_check_is_NonlinearGenomicModel_is_concrete():
    generic_assert_concrete_function(check_is_NonlinearGenomicModel)

def test_cond_check_is_NonlinearGenomicModel_is_concrete():
    generic_assert_concrete_function(cond_check_is_NonlinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_NonlinearGenomicModel(gmod):
    assert is_NonlinearGenomicModel(gmod)

def test_check_is_NonlinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_NonlinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_NonlinearGenomicModel(None, "gmod")
