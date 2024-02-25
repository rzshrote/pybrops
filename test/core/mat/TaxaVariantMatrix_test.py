import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.core.mat.TaxaVariantMatrix import check_is_TaxaVariantMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyTaxaVariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(TaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TaxaVariantMatrix_is_concrete():
    assert_function_isconcrete(check_is_TaxaVariantMatrix)

def test_check_is_TaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaVariantMatrix(None, "mat")
