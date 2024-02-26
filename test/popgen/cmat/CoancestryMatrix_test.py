import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.CoancestryMatrix import check_is_CoancestryMatrix
from .common_fixtures import *

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
    assert_class_documentation(CoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mat_asformat_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "mat_asformat")

def test_coancestry_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "coancestry")

def test_kinship_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "kinship")

def test_is_positive_semidefinite_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "is_positive_semidefinite")

def test_apply_jitter_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "apply_jitter")

def test_max_inbreeding_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "max_inbreeding")

def test_min_inbreeding_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "min_inbreeding")

def test_inverse_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "inverse")

def test_max_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "max")

def test_mean_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "mean")

def test_min_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "min")

def test_to_pandas_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "to_pandas")

def test_to_csv_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "to_csv")

def test_to_hdf5_is_abstract():
    assert_method_isabstract(CoancestryMatrix, "to_hdf5")

def test_from_pandas_is_abstract():
    assert_classmethod_isabstract(CoancestryMatrix, "from_pandas")

def test_from_csv_is_abstract():
    assert_classmethod_isabstract(CoancestryMatrix, "from_csv")

def test_from_hdf5_is_abstract():
    assert_classmethod_isabstract(CoancestryMatrix, "from_hdf5")

def test_from_gmat_is_abstract():
    assert_classmethod_isabstract(CoancestryMatrix, "from_gmat")

################################################################################
################## Test for concrete class utility functions ###################
################################################################################
def test_check_is_CoancestryMatrix_is_concrete():
    assert_function_isconcrete(check_is_CoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_CoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_CoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_CoancestryMatrix(None, "cmat")
