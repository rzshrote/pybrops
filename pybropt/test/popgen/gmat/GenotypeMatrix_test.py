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

@pytest.fixture
def vmethods(mat):
    yield [m for m in dir(mat) if m.startswith('__') is False]

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_abstract_methods(mat, vmethods):
    generic_test_abstract_methods(mat, vmethods)

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
