import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import check_is_AdditiveGeneticVarianceMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyAdditiveGeneticVarianceMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveGeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_classmethod_isabstract(AdditiveGeneticVarianceMatrix, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGeneticVarianceMatrix_is_concrete():
    assert_function_isconcrete(check_is_AdditiveGeneticVarianceMatrix)

def test_check_is_AdditiveGeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(None, "mat")
