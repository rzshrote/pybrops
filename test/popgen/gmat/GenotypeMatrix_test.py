import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyGenotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(GenotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ploidy_is_abstract():
    assert_property_isabstract(GenotypeMatrix, "ploidy")

def test_nphase_is_abstract():
    assert_property_isabstract(GenotypeMatrix, "nphase")

def test_mat_format_is_abstract():
    assert_property_isabstract(GenotypeMatrix, "mat_format")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mat_asformat_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "mat_asformat")

def test_tacount_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "tacount")

def test_tafreq_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "tafreq")

def test_acount_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "acount")

def test_afreq_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "afreq")

def test_apoly_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "apoly")

def test_maf_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "maf")

def test_meh_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "meh")

def test_gtcount_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "gtcount")

def test_gtfreq_is_abstract():
    assert_method_isabstract(GenotypeMatrix, "gtfreq")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GenotypeMatrix_is_concrete():
    assert_function_isconcrete(check_is_GenotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_GenotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GenotypeMatrix(None, "mat")
