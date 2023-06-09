import inspect
import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix
from pybrops.popgen.gmat.HaplotypeMatrix import check_is_HaplotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield HaplotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(HaplotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(HaplotypeMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ploidy_is_abstract():
    assert_abstract_property(HaplotypeMatrix, "ploidy")

def test_nphase_is_abstract():
    assert_abstract_property(HaplotypeMatrix, "nphase")

def test_mat_format_is_abstract():
    assert_abstract_property(HaplotypeMatrix, "mat_format")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_thcount_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "thcount")

def test_thfreq_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "thfreq")

def test_hcount_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "hcount")

def test_hfreq_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "hfreq")

def test_mhf_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "mhf")

def test_meh_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "meh")

def test_gtcount_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "gtcount")

def test_gtfreq_is_abstract():
    assert_abstract_method(HaplotypeMatrix, "gtfreq")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_HaplotypeMatrix_is_concrete():
    assert_concrete_function(check_is_HaplotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_HaplotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_HaplotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_HaplotypeMatrix(None, "mat")
