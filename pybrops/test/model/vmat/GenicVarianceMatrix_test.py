import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.model.vmat.GenicVarianceMatrix import check_is_GenicVarianceMatrix
from pybrops.test.model.vmat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyGenicVarianceMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GenicVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GenicVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract():
    assert_abstract_method(GenicVarianceMatrix, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenicVarianceMatrix_is_concrete():
    assert_concrete_function(check_is_GenicVarianceMatrix)

def test_check_is_GenicVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_GenicVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GenicVarianceMatrix(None, "mat")
