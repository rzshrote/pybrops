import pytest

from pybrops.opt.algo.RealOptimizationAlgorithm import RealOptimizationAlgorithm, check_is_RealOptimizationAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.test.assert_python import assert_class_isabstract, assert_method_isabstract, assert_class_documentation, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

class DummyRealOptimizationAlgorithm(RealOptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummyRealOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_RealOptimizationAlgorithm_is_semiabstract():
    assert_class_isabstract(RealOptimizationAlgorithm)

############################## Test class docstring ############################
def test_RealOptimizationAlgorithm_docstring():
    assert_class_documentation(RealOptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_RealOptimizationAlgorithm_minimize_is_abstract():
    assert_method_isabstract(RealOptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_RealOptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_RealOptimizationAlgorithm(algo, "algo")

def test_check_is_RealOptimizationAlgorithm_TypeError(algo):
    with pytest.raises(TypeError):
        check_is_RealOptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_RealOptimizationAlgorithm(DummyOptimizationAlgorithm(), "algo")
    with pytest.raises(TypeError):
        check_is_RealOptimizationAlgorithm(None, "algo")
