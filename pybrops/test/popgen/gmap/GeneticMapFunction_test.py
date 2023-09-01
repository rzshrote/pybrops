import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import check_is_GeneticMapFunction
from pybrops.test.popgen.gmap.common_fixtures import *

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
    assert_docstring(GeneticMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GeneticMapFunction, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mapfn_is_abstract():
    assert_abstract_method(GeneticMapFunction, "mapfn")

def test_invmapfn_is_abstract():
    assert_abstract_method(GeneticMapFunction, "invmapfn")

def test_rprob1g_is_abstract():
    assert_abstract_method(GeneticMapFunction, "rprob1g")

def test_rprob2g_is_abstract():
    assert_abstract_method(GeneticMapFunction, "rprob2g")

def test_rprob1p_is_abstract():
    assert_abstract_method(GeneticMapFunction, "rprob1p")

def test_rprob2p_is_abstract():
    assert_abstract_method(GeneticMapFunction, "rprob2p")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GeneticMapFunction_is_concrete():
    assert_concrete_function(check_is_GeneticMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GeneticMapFunction(gmapfn):
    with not_raises(TypeError):
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
    with pytest.raises(TypeError):
        check_is_GeneticMapFunction(None, "gmapfn")
