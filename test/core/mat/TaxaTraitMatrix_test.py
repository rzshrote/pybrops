import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import check_is_TaxaTraitMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyTaxaTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(TaxaTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TaxaTraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_TaxaTraitMatrix)

def test_check_is_TaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaTraitMatrix(None, "mat")
