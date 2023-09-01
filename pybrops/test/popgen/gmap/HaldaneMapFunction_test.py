import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction
from pybrops.popgen.gmap.HaldaneMapFunction import check_is_HaldaneMapFunction

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
    assert_docstring(HaldaneMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "__init__")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "mapfn")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "invmapfn")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "rprob1g")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "rprob2g")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "rprob1p")

def test_mapfn_is_concrete():
    assert_concrete_method(HaldaneMapFunction, "rprob2p")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

### mapfn
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

## invmapfn
def test_invmapfn_rand(mapfn):
    dist = numpy.random.uniform(0,3,100)
    recomb = mapfn.mapfn(dist)
    out = mapfn.invmapfn(recomb)
    assert numpy.allclose(out, dist, equal_nan = True)

def test_invmapfn_zero(mapfn):
    dist = numpy.repeat(0.0, 100)
    recomb = mapfn.mapfn(dist)
    out = mapfn.invmapfn(recomb)
    assert numpy.allclose(out, dist, equal_nan = True)

def test_invmapfn_inf(mapfn):
    dist = numpy.repeat(numpy.inf, 100)
    recomb = mapfn.mapfn(dist)
    out = mapfn.invmapfn(recomb)
    assert numpy.allclose(out, dist, equal_nan = True)

def test_invmapfn_NaN(mapfn):
    dist = numpy.repeat(numpy.nan, 100)
    recomb = mapfn.mapfn(dist)
    out = mapfn.invmapfn(recomb)
    assert numpy.allclose(out, dist, equal_nan = True)

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_HaldaneMapFunction_is_concrete():
    assert_concrete_function(check_is_HaldaneMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_HaldaneMapFunction(mapfn):
    with not_raises(TypeError):
        check_is_HaldaneMapFunction(mapfn, "mapfn")
    with pytest.raises(TypeError):
        check_is_HaldaneMapFunction(None, "mapfn")
