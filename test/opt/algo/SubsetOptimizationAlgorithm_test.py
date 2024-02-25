import pytest

from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm, check_is_SubsetOptimizationAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.test.assert_python import assert_class_isabstract, assert_method_isabstract, assert_class_documentation, not_raises

class DummyOptimizationAlgorithm(OptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

class DummySubsetOptimizationAlgorithm(SubsetOptimizationAlgorithm):
    def minimize(self, prob, miscout, **kwargs):
        return super().minimize(prob, miscout, **kwargs)

@pytest.fixture
def algo():
    yield DummySubsetOptimizationAlgorithm()

################### Test class abstract/concrete properties ####################
def test_SubsetOptimizationAlgorithm_is_semiabstract():
    assert_class_isabstract(SubsetOptimizationAlgorithm)

############################## Test class docstring ############################
def test_SubsetOptimizationAlgorithm_docstring():
    assert_class_documentation(SubsetOptimizationAlgorithm)

############################# Test class properties ############################

############################## Test class methods ##############################
def test_SubsetOptimizationAlgorithm_minimize_is_abstract():
    assert_method_isabstract(SubsetOptimizationAlgorithm, "minimize")

############################# Test class utilities #############################
def test_check_is_SubsetOptimizationAlgorithm(algo):
    with not_raises(Exception):
        check_is_SubsetOptimizationAlgorithm(algo, "algo")

def test_check_is_SubsetOptimizationAlgorithm_TypeError(algo):
    with pytest.raises(TypeError):
        check_is_SubsetOptimizationAlgorithm(object(), "algo")
    with pytest.raises(TypeError):
        check_is_SubsetOptimizationAlgorithm(DummyOptimizationAlgorithm(), "algo")
    with pytest.raises(TypeError):
        check_is_SubsetOptimizationAlgorithm(None, "algo")
