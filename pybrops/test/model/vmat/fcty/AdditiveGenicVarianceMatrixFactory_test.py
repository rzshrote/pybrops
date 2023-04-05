import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.model.vmat.fcty.AdditiveGenicVarianceMatrixFactory import AdditiveGenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.AdditiveGenicVarianceMatrixFactory import check_is_AdditiveGenicVarianceMatrixFactory

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield AdditiveGenicVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(AdditiveGenicVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveGenicVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract(fcty):
    assert_abstract_method(fcty, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGenicVarianceMatrixFactory_is_concrete():
    assert_concrete_function(check_is_AdditiveGenicVarianceMatrixFactory)

def test_check_is_AdditiveGenicVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_AdditiveGenicVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_AdditiveGenicVarianceMatrixFactory(None, "fcty")
