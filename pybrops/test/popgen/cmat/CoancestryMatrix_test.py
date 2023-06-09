import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.CoancestryMatrix import check_is_CoancestryMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def cmat():
    yield CoancestryMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(CoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(CoancestryMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_coancestry_is_abstract():
    assert_abstract_method(CoancestryMatrix, "coancestry")

def test_kinship_is_abstract():
    assert_abstract_method(CoancestryMatrix, "kinship")

def test_is_positive_semidefinite_is_abstract():
    assert_abstract_method(CoancestryMatrix, "is_positive_semidefinite")

def test_apply_jitter_is_abstract():
    assert_abstract_method(CoancestryMatrix, "apply_jitter")

def test_from_gmat_is_abstract():
    assert_abstract_method(CoancestryMatrix, "from_gmat")

################################################################################
################## Test for concrete class utility functions ###################
################################################################################
def test_check_is_CoancestryMatrix_is_concrete():
    assert_concrete_function(check_is_CoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_CoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_CoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_CoancestryMatrix(None, "cmat")
