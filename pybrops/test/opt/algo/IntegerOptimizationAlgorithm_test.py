import pytest

from pybrops.opt.algo.IntegerOptimizationAlgorithm import IntegerOptimizationAlgorithm, check_is_IntegerOptimizationAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.test.assert_python import assert_abstract_class, assert_abstract_method, assert_docstring, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

class DummyIntegerOptimizationAlgorithm(IntegerOptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummyIntegerOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_IntegerOptimizationAlgorithm_is_semiabstract():
    assert_abstract_class(IntegerOptimizationAlgorithm)

############################## Test class docstring ############################
def test_IntegerOptimizationAlgorithm_docstring():
    assert_docstring(IntegerOptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_IntegerOptimizationAlgorithm_minimize_is_abstract():
    assert_abstract_method(IntegerOptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_IntegerOptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_IntegerOptimizationAlgorithm(algo, "algo")

def test_check_is_IntegerOptimizationAlgorithm_TypeError(algo):
    with pytest.raises(TypeError):
        check_is_IntegerOptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_IntegerOptimizationAlgorithm(DummyOptimizationAlgorithm(), "algo")
    with pytest.raises(TypeError):
        check_is_IntegerOptimizationAlgorithm(None, "algo")