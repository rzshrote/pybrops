import os
from matplotlib import pyplot
import pytest
from pybrops.opt.algo.NSGA2BinaryGeneticAlgorithm import NSGA2BinaryGeneticAlgorithm
from pybrops.opt.soln.BinarySolution import BinarySolution
from pybrops.test.assert_python import assert_class_isconcrete, assert_method_isconcrete, assert_property_isconcrete, assert_class_documentation, not_raises
from .common_fixtures import *

@pytest.fixture
def algo(
        common_ngen,
        common_pop_size,
        common_rng
    ):
    out = NSGA2BinaryGeneticAlgorithm(
        ngen=common_ngen,
        pop_size=common_pop_size,
        rng=common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_NSGA2BinaryGeneticAlgorithm_is_concrete():
    assert_class_isconcrete(NSGA2BinaryGeneticAlgorithm)

############################## Test class docstring ############################
def test_NSGA2BinaryGeneticAlgorithm_docstring():
    assert_class_documentation(NSGA2BinaryGeneticAlgorithm)

############################# Test class properties ############################

### ngen ###
def test_NSGA2BinaryGeneticAlgorithm_ngen_is_concrete():
    assert_property_isconcrete(NSGA2BinaryGeneticAlgorithm, "ngen")

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
def test_NSGA2BinaryGeneticAlgorithm_pop_size_is_concrete():
    assert_property_isconcrete(NSGA2BinaryGeneticAlgorithm, "pop_size")

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
def test_NSGA2BinaryGeneticAlgorithm_rng_is_concrete():
    assert_property_isconcrete(NSGA2BinaryGeneticAlgorithm, "rng")

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
def test_NSGA2BinaryGeneticAlgorithm_minimize_is_concrete():
    assert_method_isconcrete(NSGA2BinaryGeneticAlgorithm, "minimize")

def test_minimize_single(algo, common_sum_prob_binary_single):
    with pytest.raises(TypeError):
        algo.minimize(common_sum_prob_binary_single)

def test_minimize_multi(algo, common_sum_prob_binary_multi):
    soln = algo.minimize(common_sum_prob_binary_multi)
    assert isinstance(soln, BinarySolution)

    # make a directory if needed
    if not os.path.isdir("optsoln"):
        os.makedirs("optsoln")

    # save solutions
    numpy.savetxt(
        "optsoln/NSGA2BinaryGeneticAlgorithm_sum_prob_multi_soln_decn.test.txt",
        soln.soln_decn,
        fmt="%i"
    )

    # make graph
    fig = pyplot.figure()
    ax = pyplot.axes()
    ax.scatter(
        soln.soln_obj[:,0], 
        soln.soln_obj[:,1]
    )
    ax.set_xlabel("Objective 1")
    ax.set_ylabel("Objective 2")
    ax.set_title("NSGA2 Binary Genetic Algorithm Test Pareto Frontier")
    pyplot.savefig("optsoln/NSGA2BinaryGeneticAlgorithm_sum_prob_multi_soln_decn.test.png", dpi = 250)
    pyplot.close()

############################# Test class utilities #############################
