import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import is_NonlinearGenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import check_is_NonlinearGenomicModel

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
    assert_docstring(NonlinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(NonlinearGenomicModel, "__init__")

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
    assert_concrete_function(is_NonlinearGenomicModel)

def test_check_is_NonlinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_NonlinearGenomicModel)

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
