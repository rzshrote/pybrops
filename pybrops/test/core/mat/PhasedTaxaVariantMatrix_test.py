import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.core.mat.PhasedTaxaVariantMatrix import is_PhasedTaxaVariantMatrix
from pybrops.core.mat.PhasedTaxaVariantMatrix import check_is_PhasedTaxaVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PhasedTaxaVariantMatrix()

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
def test_is_PhasedTaxaVariantMatrix_is_concrete():
    assert_concrete_function(is_PhasedTaxaVariantMatrix)

def test_is_PhasedTaxaVariantMatrix(mat):
    assert is_PhasedTaxaVariantMatrix(mat)

def test_check_is_PhasedTaxaVariantMatrix_is_concrete():
    assert_concrete_function(check_is_PhasedTaxaVariantMatrix)

def test_check_is_PhasedTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(None, "mat")
