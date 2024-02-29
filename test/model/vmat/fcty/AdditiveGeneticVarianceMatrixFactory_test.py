import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import AdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import check_is_AdditiveGeneticVarianceMatrixFactory
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield DummyAdditiveGeneticVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveGeneticVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_method_isabstract(AdditiveGeneticVarianceMatrixFactory, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGeneticVarianceMatrixFactory_is_concrete():
    assert_function_isconcrete(check_is_AdditiveGeneticVarianceMatrixFactory)

def test_check_is_AdditiveGeneticVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrixFactory(None, "fcty")
