import os
import pytest
from pybrops.opt.algo.BinaryGeneticAlgorithm import BinaryGeneticAlgorithm
from pybrops.opt.soln.BinarySolution import BinarySolution
from pybrops.test.assert_python import assert_concrete_class, assert_concrete_method, assert_concrete_property, assert_docstring, not_raises
from .common_fixtures import *

@pytest.fixture
def algo(
        common_ngen,
        common_pop_size,
        common_rng
    ):
    out = BinaryGeneticAlgorithm(
        ngen=common_ngen,
        pop_size=common_pop_size,
        rng=common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_BinaryGeneticAlgorithm_is_concrete():
    assert_concrete_class(BinaryGeneticAlgorithm)

############################## Test class docstring ############################
def test_BinaryGeneticAlgorithm_docstring():
    assert_docstring(BinaryGeneticAlgorithm)

############################# Test class properties ############################

### ngen ###
def test_BinaryGeneticAlgorithm_ngen_is_concrete():
    assert_concrete_property(BinaryGeneticAlgorithm, "ngen")

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
def test_BinaryGeneticAlgorithm_pop_size_is_concrete():
    assert_concrete_property(BinaryGeneticAlgorithm, "pop_size")

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
def test_BinaryGeneticAlgorithm_rng_is_concrete():
    assert_concrete_property(BinaryGeneticAlgorithm, "rng")

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
def test_BinaryGeneticAlgorithm_minimize_is_concrete():
    assert_concrete_method(BinaryGeneticAlgorithm, "minimize")

def test_minimize_single(algo, common_sum_prob_binary_single):
    soln = algo.minimize(common_sum_prob_binary_single)
    assert isinstance(soln, BinarySolution)

    # make a directory if needed
    if not os.path.isdir("optsoln"):
        os.makedirs("optsoln")

    # save solutions
    numpy.savetxt(
        "optsoln/BinaryGeneticAlgorithm_sum_prob_single_soln_decn.test.txt",
        soln.soln_decn,
        fmt="%i"
    )

def test_minimize_multi(algo, common_sum_prob_binary_multi):
    with pytest.raises(TypeError):
        algo.minimize(common_sum_prob_binary_multi)

############################# Test class utilities #############################