import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import check_is_GeneticVarianceMatrix
from .common_fixtures import *

################################ Test fixtures #################################
@pytest.fixture
def mat():
    yield DummyGeneticVarianceMatrix()

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(GeneticVarianceMatrix)

############################# Test concrete methods ############################

########################### Test abstract properties ###########################

### epgc
def test_epgc_is_abstract():
    assert_property_isabstract(GeneticVarianceMatrix, "epgc")

############################# Test abstract methods ############################

### from_gmod
def test_from_gmod_is_abstract():
    assert_classmethod_isabstract(GeneticVarianceMatrix, "from_gmod")

######################### Test class utility functions #########################
def test_check_is_GeneticVarianceMatrix_is_concrete():
    assert_function_isconcrete(check_is_GeneticVarianceMatrix)

def test_check_is_GeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_GeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GeneticVarianceMatrix(None, "mat")
