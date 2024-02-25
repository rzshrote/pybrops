import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import check_is_GenicVarianceMatrixFactory
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield DummyGenicVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(GenicVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(GenicVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_method_isabstract(GenicVarianceMatrixFactory, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenicVarianceMatrixFactory_is_concrete():
    assert_function_isconcrete(check_is_GenicVarianceMatrixFactory)

def test_check_is_GenicVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_GenicVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_GenicVarianceMatrixFactory(None, "fcty")
