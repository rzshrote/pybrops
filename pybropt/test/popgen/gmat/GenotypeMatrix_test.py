import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.gmat import GenotypeMatrix
from pybropt.popgen.gmat import is_GenotypeMatrix
from pybropt.popgen.gmat import check_is_GenotypeMatrix
from pybropt.popgen.gmat import cond_check_is_GenotypeMatrix


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield GenotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GenotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GenotypeMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ploidy_is_abstract():
    generic_assert_abstract_property(GenotypeMatrix, "ploidy")

def test_nphase_is_abstract():
    generic_assert_abstract_property(GenotypeMatrix, "nphase")

def test_mat_format_is_abstract():
    generic_assert_abstract_property(GenotypeMatrix, "mat_format")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mat_asformat_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "mat_asformat")

def test_tacount_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "tacount")

def test_tafreq_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "tafreq")

def test_acount_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "acount")

def test_afreq_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "afreq")

def test_maf_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "maf")

def test_mehe_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "mehe")

def test_gtcount_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "gtcount")

def test_gtfreq_is_abstract():
    generic_assert_abstract_method(GenotypeMatrix, "gtfreq")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_GenotypeMatrix_is_concrete():
    generic_assert_concrete_function(is_GenotypeMatrix)

def test_check_is_GenotypeMatrix_is_concrete():
    generic_assert_concrete_function(check_is_GenotypeMatrix)

def test_cond_check_is_GenotypeMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_GenotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GenotypeMatrix(mat):
    assert is_GenotypeMatrix(mat)

def test_check_is_GenotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_GenotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GenotypeMatrix(None, "mat")
