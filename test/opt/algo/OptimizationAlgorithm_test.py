import pytest

from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.test.assert_python import assert_class_isabstract, assert_method_isabstract, assert_class_documentation, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummyOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_OptimizationAlgorithm_is_semiabstract():
    assert_class_isabstract(OptimizationAlgorithm)

############################## Test class docstring ############################
def test_OptimizationAlgorithm_docstring():
    assert_class_documentation(OptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_OptimizationAlgorithm_minimize_is_abstract():
    assert_method_isabstract(OptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_OptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_OptimizationAlgorithm(algo, "algo")
    with pytest.raises(TypeError):
        check_is_OptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_OptimizationAlgorithm(None, "algo")
