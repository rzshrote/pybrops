import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import check_is_GeneticVarianceMatrixFactory
from .common_fixtures import *

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
    assert_class_documentation(GeneticVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(GeneticVarianceMatrixFactory, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_gmod_is_abstract():
    assert_method_isabstract(GeneticVarianceMatrixFactory, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticVarianceMatrixFactory_is_concrete():
    assert_function_isconcrete(check_is_GeneticVarianceMatrixFactory)

def test_check_is_GeneticVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_GeneticVarianceMatrixFactory(None, "fcty")
