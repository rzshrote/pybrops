import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import is_GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import check_is_GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import cond_check_is_GeneticMapFunction


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmapfn():
    yield GeneticMapFunction()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GeneticMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GeneticMapFunction, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mapfn_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "mapfn")

def test_invmapfn_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "invmapfn")

def test_rprob1g_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "rprob1g")

def test_rprob2g_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "rprob2g")

def test_rprob1p_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "rprob1p")

def test_rprob2p_is_abstract():
    generic_assert_abstract_method(GeneticMapFunction, "rprob2p")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_GeneticMapFunction_is_concrete():
    generic_assert_concrete_function(is_GeneticMapFunction)

def test_check_is_GeneticMapFunction_is_concrete():
    generic_assert_concrete_function(check_is_GeneticMapFunction)

def test_cond_check_is_GeneticMapFunction_is_concrete():
    generic_assert_concrete_function(cond_check_is_GeneticMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GeneticMapFunction(gmapfn):
    assert is_GeneticMapFunction(gmapfn)

def test_check_is_GeneticMapFunction(gmapfn):
    with not_raises(TypeError):
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
    with pytest.raises(TypeError):
        check_is_GeneticMapFunction(None, "gmapfn")