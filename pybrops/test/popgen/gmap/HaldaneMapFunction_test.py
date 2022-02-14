import pytest
import numpy

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.HaldaneMapFunction import is_HaldaneMapFunction
from pybrops.popgen.gmap.HaldaneMapFunction import check_is_HaldaneMapFunction
from pybrops.popgen.gmap.HaldaneMapFunction import cond_check_is_HaldaneMapFunction

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmap(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")

@pytest.fixture
def mapfn():
    yield HaldaneMapFunction()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(HaldaneMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(HaldaneMapFunction, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_mapfn_rand(mapfn):
    d = numpy.random.uniform(0,3,100)
    r = 0.5 * (1.0 - numpy.exp(-2.0 * d))
    b = mapfn.mapfn(d)
    assert numpy.all(b == r)

def test_mapfn_zero(mapfn):
    d = numpy.repeat(0.0, 100)
    r = numpy.repeat(0.0, 100)
    b = mapfn.mapfn(d)
    assert numpy.all(b == r)

def test_mapfn_inf(mapfn):
    d = numpy.repeat(numpy.inf, 100)
    r = numpy.repeat(0.5, 100)
    b = mapfn.mapfn(d)
    assert numpy.all(b == r)

def test_mapfn_NaN(mapfn):
    d = numpy.repeat(numpy.nan, 100)
    r = numpy.repeat(numpy.nan, 100)
    b = mapfn.mapfn(d)
    assert numpy.all(numpy.isnan(b) == numpy.isnan(r))

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_HaldaneMapFunction_is_concrete():
    generic_assert_concrete_function(is_HaldaneMapFunction)

def test_check_is_HaldaneMapFunction_is_concrete():
    generic_assert_concrete_function(check_is_HaldaneMapFunction)

def test_cond_check_is_HaldaneMapFunction_is_concrete():
    generic_assert_concrete_function(cond_check_is_HaldaneMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_HaldaneMapFunction(mapfn):
    assert is_HaldaneMapFunction(mapfn)

def test_check_is_HaldaneMapFunction(mapfn):
    with not_raises(TypeError):
        check_is_HaldaneMapFunction(mapfn, "mapfn")
    with pytest.raises(TypeError):
        check_is_HaldaneMapFunction(None, "mapfn")
