import os
import pytest
from pybrops.opt.algo.RealGeneticAlgorithm import RealGeneticAlgorithm
from pybrops.opt.soln.RealSolution import RealSolution
from pybrops.test.assert_python import assert_concrete_class, assert_concrete_method, assert_concrete_property, assert_docstring, not_raises
from pybrops.test.opt.algo.common_fixtures import *

@pytest.fixture
def algo(
        common_ngen,
        common_pop_size,
        common_rng
    ):
    out = RealGeneticAlgorithm(
        ngen=common_ngen,
        pop_size=common_pop_size,
        rng=common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_RealGeneticAlgorithm_is_concrete():
    assert_concrete_class(RealGeneticAlgorithm)

############################## Test class docstring ############################
def test_RealGeneticAlgorithm_docstring():
    assert_docstring(RealGeneticAlgorithm)

############################# Test class properties ############################

### ngen ###
def test_RealGeneticAlgorithm_ngen_is_concrete():
    assert_concrete_property(RealGeneticAlgorithm, "ngen")

def test_ngen_fget(algo, common_ngen):
    assert algo.ngen == common_ngen

def test_ngen_fset(algo, common_ngen):
    with not_raises(Exception):
        algo.ngen = int(common_ngen)
    with not_raises(Exception):
        algo.ngen = numpy.int8(1)
    with not_raises(Exception):
        algo.ngen = numpy.int16(1)
    with not_raises(Exception):
        algo.ngen = numpy.int32(1)
    with not_raises(Exception):
        algo.ngen = numpy.int64(1)

def test_ngen_fset_TypeError(algo, common_ngen):
    with pytest.raises(TypeError):
        algo.ngen = object()
    with pytest.raises(TypeError):
        algo.ngen = None
    with pytest.raises(TypeError):
        algo.ngen = str(common_ngen)
    with pytest.raises(TypeError):
        algo.ngen = numpy.array([common_ngen], dtype=int)

def test_ngen_fset_ValueError(algo, common_ngen):
    with pytest.raises(ValueError):
        algo.ngen = -common_ngen
    with pytest.raises(ValueError):
        algo.ngen = int(-1)
    with pytest.raises(ValueError):
        algo.ngen = int(0)

def test_ngen_fdel(algo):
    with pytest.raises(AttributeError):
        del algo.ngen

### pop_size ###
def test_RealGeneticAlgorithm_pop_size_is_concrete():
    assert_concrete_property(RealGeneticAlgorithm, "pop_size")

def test_pop_size_fget(algo, common_pop_size):
    assert algo.pop_size == common_pop_size

def test_pop_size_fset(algo, common_pop_size):
    with not_raises(Exception):
        algo.pop_size = int(common_pop_size)
    with not_raises(Exception):
        algo.pop_size = numpy.int8(1)
    with not_raises(Exception):
        algo.pop_size = numpy.int16(1)
    with not_raises(Exception):
        algo.pop_size = numpy.int32(1)
    with not_raises(Exception):
        algo.pop_size = numpy.int64(1)

def test_pop_size_fset_TypeError(algo, common_pop_size):
    with pytest.raises(TypeError):
        algo.pop_size = object()
    with pytest.raises(TypeError):
        algo.pop_size = None
    with pytest.raises(TypeError):
        algo.pop_size = str(common_pop_size)
    with pytest.raises(TypeError):
        algo.pop_size = numpy.array([common_pop_size], dtype=int)

def test_pop_size_fset_ValueError(algo, common_pop_size):
    with pytest.raises(ValueError):
        algo.pop_size = -common_pop_size
    with pytest.raises(ValueError):
        algo.pop_size = int(-1)
    with pytest.raises(ValueError):
        algo.pop_size = int(0)

def test_pop_size_fdel(algo):
    with pytest.raises(AttributeError):
        del algo.pop_size

### rng ###
def test_RealGeneticAlgorithm_rng_is_concrete():
    assert_concrete_property(RealGeneticAlgorithm, "rng")

def test_rng_fget(algo, common_rng):
    assert algo.rng == common_rng

def test_rng_fset(algo, common_rng):
    with not_raises(Exception):
        algo.rng = common_rng
    with not_raises(Exception):
        algo.rng = None

def test_rng_fset_TypeError(algo):
    with pytest.raises(TypeError):
        algo.rng = object()
    with pytest.raises(TypeError):
        algo.rng = "string"
    with pytest.raises(TypeError):
        algo.rng = numpy.array([], dtype=int)

def test_rng_fdel(algo):
    with pytest.raises(AttributeError):
        del algo.rng

############################## Test class methods ##############################

### minimize ###
def test_RealGeneticAlgorithm_minimize_is_concrete():
    assert_concrete_method(RealGeneticAlgorithm, "minimize")

def test_minimize_single(algo, common_sum_prob_real_single):
    soln = algo.minimize(common_sum_prob_real_single)
    assert isinstance(soln, RealSolution)

    # make a directory if needed
    if not os.path.isdir("optsoln"):
        os.makedirs("optsoln")

    # save solutions
    numpy.savetxt(
        "optsoln/RealGeneticAlgorithm_sum_prob_single_soln_decn.test.txt",
        soln.soln_decn,
        fmt="%.5e"
    )

def test_minimize_multi(algo, common_sum_prob_real_multi):
    with pytest.raises(TypeError):
        algo.minimize(common_sum_prob_real_multi)

############################# Test class utilities #############################
