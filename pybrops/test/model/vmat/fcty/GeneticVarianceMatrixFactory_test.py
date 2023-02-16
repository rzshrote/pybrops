import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import check_is_GeneticVarianceMatrixFactory

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield GeneticVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GeneticVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GeneticVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract(fcty):
    generic_assert_abstract_method(fcty, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticVarianceMatrixFactory_is_concrete():
    generic_assert_concrete_function(check_is_GeneticVarianceMatrixFactory)

def test_check_is_GeneticVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(None, "fcty")
