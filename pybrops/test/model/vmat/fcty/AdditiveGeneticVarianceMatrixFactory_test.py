import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import AdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import check_is_AdditiveGeneticVarianceMatrixFactory

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield AdditiveGeneticVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(AdditiveGeneticVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveGeneticVarianceMatrixFactory, "__init__")

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
def test_check_is_AdditiveGeneticVarianceMatrixFactory_is_concrete():
    assert_concrete_function(check_is_AdditiveGeneticVarianceMatrixFactory)

def test_check_is_AdditiveGeneticVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrixFactory(None, "fcty")