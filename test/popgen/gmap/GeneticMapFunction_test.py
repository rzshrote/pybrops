import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import check_is_GeneticMapFunction
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmapfn():
    yield DummyGeneticMapFunction()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(GeneticMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mapfn_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "mapfn")

def test_invmapfn_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "invmapfn")

def test_rprob1g_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "rprob1g")

def test_rprob2g_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "rprob2g")

def test_rprob1p_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "rprob1p")

def test_rprob2p_is_abstract():
    assert_method_isabstract(GeneticMapFunction, "rprob2p")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GeneticMapFunction_is_concrete():
    assert_function_isconcrete(check_is_GeneticMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticMapFunction(gmapfn):
    with not_raises(TypeError):
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
    with pytest.raises(TypeError):
        check_is_GeneticMapFunction(None, "gmapfn")
