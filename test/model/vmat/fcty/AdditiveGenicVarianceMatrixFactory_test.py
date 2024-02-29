import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.fcty.AdditiveGenicVarianceMatrixFactory import AdditiveGenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.AdditiveGenicVarianceMatrixFactory import check_is_AdditiveGenicVarianceMatrixFactory
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield DummyAdditiveGenicVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveGenicVarianceMatrixFactory)

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
    assert_method_isabstract(AdditiveGenicVarianceMatrixFactory, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGenicVarianceMatrixFactory_is_concrete():
    assert_function_isconcrete(check_is_AdditiveGenicVarianceMatrixFactory)

def test_check_is_AdditiveGenicVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_AdditiveGenicVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_AdditiveGenicVarianceMatrixFactory(None, "fcty")
