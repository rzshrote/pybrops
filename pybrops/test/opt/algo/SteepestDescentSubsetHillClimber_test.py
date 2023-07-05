import os
import pytest
from pybrops.opt.algo.SteepestDescentSubsetHillClimber import SteepestDescentSubsetHillClimber
from pybrops.opt.soln.SubsetSolution import SubsetSolution
from pybrops.test.assert_python import assert_concrete_class, assert_concrete_method, assert_concrete_property, assert_docstring, not_raises
from pybrops.test.opt.algo.common_fixtures import *

@pytest.fixture
def algo(
        common_rng
    ):
    out = SteepestDescentSubsetHillClimber(
        rng=common_rng
    )
    yield out

################### Test class abstract/concrete properties ####################
def test_SteepestDescentSubsetHillClimber_is_concrete():
    assert_concrete_class(SteepestDescentSubsetHillClimber)

############################## Test class docstring ############################
def test_SteepestDescentSubsetHillClimber_docstring():
    assert_docstring(SteepestDescentSubsetHillClimber)

############################# Test class properties ############################

### rng ###
def test_SteepestDescentSubsetHillClimber_rng_is_concrete():
    assert_concrete_property(SteepestDescentSubsetHillClimber, "rng")

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
def test_SteepestDescentSubsetHillClimber_minimize_is_concrete():
    assert_concrete_method(SteepestDescentSubsetHillClimber, "minimize")

def test_minimize_single(algo, common_sum_prob_subset_single):
    soln = algo.minimize(common_sum_prob_subset_single)
    assert isinstance(soln, SubsetSolution)

    # make a directory if needed
    if not os.path.isdir("optsoln"):
        os.makedirs("optsoln")

    # save solutions
    numpy.savetxt(
        "optsoln/SteepestDescentSubsetHillClimber_sum_prob_single_soln_decn.test.txt",
        soln.soln_decn,
        fmt="%i"
    )

def test_minimize_multi(algo, common_sum_prob_subset_multi):
    with pytest.raises(TypeError):
        algo.minimize(common_sum_prob_subset_multi)

############################# Test class utilities #############################
