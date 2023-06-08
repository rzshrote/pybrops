import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import check_is_TaxaTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield TaxaTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(TaxaTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(TaxaTraitMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TaxaTraitMatrix_is_concrete():
    assert_concrete_function(check_is_TaxaTraitMatrix)

def test_check_is_TaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaTraitMatrix(None, "mat")
