import pytest

from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.test.assert_python import assert_abstract_class, assert_abstract_method, assert_docstring, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummyOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_OptimizationAlgorithm_is_semiabstract():
    assert_abstract_class(OptimizationAlgorithm)

############################## Test class docstring ############################
def test_OptimizationAlgorithm_docstring():
    assert_docstring(OptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_OptimizationAlgorithm_minimize_is_abstract():
    assert_abstract_method(OptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_OptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_OptimizationAlgorithm(algo, "algo")
    with pytest.raises(TypeError):
        check_is_OptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_OptimizationAlgorithm(None, "algo")
