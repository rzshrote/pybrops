import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import check_is_NonlinearGenomicModel
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmod():
    yield DummyNonlinearGenomicModel()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(NonlinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_NonlinearGenomicModel_is_concrete():
    assert_function_isconcrete(check_is_NonlinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_NonlinearGenomicModel(gmod):
    with not_raises(TypeError):
        check_is_NonlinearGenomicModel(gmod, "gmod")
    with pytest.raises(TypeError):
        check_is_NonlinearGenomicModel(None, "gmod")
