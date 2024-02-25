import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

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
    assert_class_documentation(PhasedTaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhasedTaxaVariantMatrix_is_concrete():
    assert_function_isconcrete(check_is_PhasedTaxaVariantMatrix)

def test_check_is_PhasedTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(None, "mat")
