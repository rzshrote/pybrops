import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import check_is_GenicVarianceMatrixFactory

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield GenicVarianceMatrixFactory()

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
def test_from_algmod_is_abstract(fcty):
    assert_abstract_method(fcty, "from_gmod")

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
