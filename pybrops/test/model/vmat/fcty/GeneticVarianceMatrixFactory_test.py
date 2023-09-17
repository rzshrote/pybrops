import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import check_is_GeneticVarianceMatrixFactory
from pybrops.test.model.vmat.fcty.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def fcty():
    yield DummyGeneticVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GeneticVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GeneticVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_gmod_is_abstract():
    assert_abstract_method(GeneticVarianceMatrixFactory, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticVarianceMatrixFactory_is_concrete():
    assert_concrete_function(check_is_GeneticVarianceMatrixFactory)

def test_check_is_GeneticVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(None, "fcty")
