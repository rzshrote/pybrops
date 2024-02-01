import pytest

from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm, check_is_BinaryOptimizationAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.test.assert_python import assert_abstract_class, assert_abstract_method, assert_docstring, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

class DummyBinaryOptimizationAlgorithm(BinaryOptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummyBinaryOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_BinaryOptimizationAlgorithm_is_semiabstract():
    assert_abstract_class(BinaryOptimizationAlgorithm)

############################## Test class docstring ############################
def test_BinaryOptimizationAlgorithm_docstring():
    assert_docstring(BinaryOptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_BinaryOptimizationAlgorithm_minimize_is_abstract():
    assert_abstract_method(BinaryOptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_BinaryOptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_BinaryOptimizationAlgorithm(algo, "algo")

def test_check_is_BinaryOptimizationAlgorithm_TypeError(algo):
    with pytest.raises(TypeError):
        check_is_BinaryOptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_BinaryOptimizationAlgorithm(DummyOptimizationAlgorithm(), "algo")
    with pytest.raises(TypeError):
        check_is_BinaryOptimizationAlgorithm(None, "algo")
