import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

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
    assert_docstring(TaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(TaxaVariantMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TaxaVariantMatrix_is_concrete():
    assert_concrete_function(check_is_TaxaVariantMatrix)

def test_check_is_TaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaVariantMatrix(None, "mat")