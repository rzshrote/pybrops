import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.CoancestryMatrix import check_is_CoancestryMatrix
from pybrops.test.popgen.cmat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def cmat():
    yield DummyCoancestryMatrix()

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
def test_mat_asformat_is_abstract():
    assert_abstract_method(CoancestryMatrix, "mat_asformat")

def test_coancestry_is_abstract():
    assert_abstract_method(CoancestryMatrix, "coancestry")

def test_kinship_is_abstract():
    assert_abstract_method(CoancestryMatrix, "kinship")

def test_is_positive_semidefinite_is_abstract():
    assert_abstract_method(CoancestryMatrix, "is_positive_semidefinite")

def test_apply_jitter_is_abstract():
    assert_abstract_method(CoancestryMatrix, "apply_jitter")

def test_max_inbreeding_is_abstract():
    assert_abstract_method(CoancestryMatrix, "max_inbreeding")

def test_min_inbreeding_is_abstract():
    assert_abstract_method(CoancestryMatrix, "min_inbreeding")

def test_inverse_is_abstract():
    assert_abstract_method(CoancestryMatrix, "inverse")

def test_max_is_abstract():
    assert_abstract_method(CoancestryMatrix, "max")

def test_mean_is_abstract():
    assert_abstract_method(CoancestryMatrix, "mean")

def test_min_is_abstract():
    assert_abstract_method(CoancestryMatrix, "min")

def test_to_pandas_is_abstract():
    assert_abstract_method(CoancestryMatrix, "to_pandas")

def test_to_csv_is_abstract():
    assert_abstract_method(CoancestryMatrix, "to_csv")

def test_to_hdf5_is_abstract():
    assert_abstract_method(CoancestryMatrix, "to_hdf5")

def test_from_pandas_is_abstract():
    assert_abstract_method(CoancestryMatrix, "from_pandas")

def test_from_csv_is_abstract():
    assert_abstract_method(CoancestryMatrix, "from_csv")

def test_from_hdf5_is_abstract():
    assert_abstract_method(CoancestryMatrix, "from_hdf5")

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
