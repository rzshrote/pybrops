import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix
from pybrops.popgen.gmat.HaplotypeMatrix import check_is_HaplotypeMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyHaplotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(HaplotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(HaplotypeMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ploidy_is_abstract():
    assert_property_isabstract(HaplotypeMatrix, "ploidy")

def test_nphase_is_abstract():
    assert_property_isabstract(HaplotypeMatrix, "nphase")

def test_mat_format_is_abstract():
    assert_property_isabstract(HaplotypeMatrix, "mat_format")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_thcount_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "thcount")

def test_thfreq_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "thfreq")

def test_hcount_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "hcount")

def test_hfreq_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "hfreq")

def test_mhf_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "mhf")

def test_meh_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "meh")

def test_gtcount_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "gtcount")

def test_gtfreq_is_abstract():
    assert_method_isabstract(HaplotypeMatrix, "gtfreq")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_HaplotypeMatrix_is_concrete():
    assert_function_isconcrete(check_is_HaplotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_HaplotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_HaplotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_HaplotypeMatrix(None, "mat")
