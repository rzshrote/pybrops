import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.model.vmat.GenicVarianceMatrix import check_is_GenicVarianceMatrix
from .common_fixtures import *

################################ Test fixtures #################################
@pytest.fixture
def mat():
    yield DummyGenicVarianceMatrix()

############################## Test class docstring ############################
def test_class_docstring():
    assert_class_documentation(GenicVarianceMatrix)

############################# Test concrete methods ############################

########################### Test abstract properties ###########################

### epgc
def test_epgc_is_abstract():
    assert_property_isabstract(GenicVarianceMatrix, "epgc")

############################# Test abstract methods ############################

### from_gmod
def test_from_gmod_is_abstract():
    assert_classmethod_isabstract(GenicVarianceMatrix, "from_gmod")

######################### Test class utility functions #########################
def test_check_is_GenicVarianceMatrix_is_concrete():
    assert_function_isconcrete(check_is_GenicVarianceMatrix)

def test_check_is_GenicVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_GenicVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GenicVarianceMatrix(None, "mat")
