import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import check_is_GenicVarianceMatrixFactory
from pybrops.test.model.vmat.fcty.common_fixtures import *

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
    assert_docstring(GenicVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GenicVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_abstract_method(GenicVarianceMatrixFactory, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenicVarianceMatrixFactory_is_concrete():
    assert_concrete_function(check_is_GenicVarianceMatrixFactory)

def test_check_is_GenicVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_GenicVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_GenicVarianceMatrixFactory(None, "fcty")
