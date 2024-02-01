import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.core.mat.PhasedTaxaVariantMatrix import check_is_PhasedTaxaVariantMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyPhasedTaxaVariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(PhasedTaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PhasedTaxaVariantMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhasedTaxaVariantMatrix_is_concrete():
    assert_concrete_function(check_is_PhasedTaxaVariantMatrix)

def test_check_is_PhasedTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(None, "mat")
