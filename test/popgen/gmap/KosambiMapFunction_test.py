import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.KosambiMapFunction import KosambiMapFunction
from pybrops.popgen.gmap.KosambiMapFunction import check_is_KosambiMapFunction

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def gmap(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")

@pytest.fixture
def mapfn():
    yield KosambiMapFunction()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(KosambiMapFunction)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "__init__")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "mapfn")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "invmapfn")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "rprob1g")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "rprob2g")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "rprob1p")

def test_mapfn_is_concrete():
    assert_method_isconcrete(KosambiMapFunction, "rprob2p")

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
    dist = numpy.random.uniform(0,3,100)
    recomb = 0.5 * numpy.tanh(2.0 * dist)
    out = mapfn.mapfn(dist)
    assert numpy.all(out == recomb)

def test_mapfn_zero(mapfn):
    dist = numpy.repeat(0.0, 100)
    recomb = numpy.repeat(0.0, 100)
    out = mapfn.mapfn(dist)
    assert numpy.all(out == recomb)

def test_mapfn_inf(mapfn):
    dist = numpy.repeat(numpy.inf, 100)
    recomb = numpy.repeat(0.5, 100)
    out = mapfn.mapfn(dist)
    assert numpy.all(out == recomb)

def test_mapfn_NaN(mapfn):
    dist = numpy.repeat(numpy.nan, 100)
    recomb = numpy.repeat(numpy.nan, 100)
    out = mapfn.mapfn(dist)
    assert numpy.all(numpy.isnan(out) == numpy.isnan(recomb))

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
def test_check_is_KosambiMapFunction_is_concrete():
    assert_function_isconcrete(check_is_KosambiMapFunction)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_KosambiMapFunction(mapfn):
    with not_raises(TypeError):
        check_is_KosambiMapFunction(mapfn, "mapfn")
    with pytest.raises(TypeError):
        check_is_KosambiMapFunction(None, "mapfn")
