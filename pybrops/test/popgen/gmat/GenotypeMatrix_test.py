import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix
from pybrops.test.popgen.gmat.common_fixtures import *

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
    assert_docstring(GenotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GenotypeMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_ploidy_is_abstract():
    assert_abstract_property(GenotypeMatrix, "ploidy")

def test_nphase_is_abstract():
    assert_abstract_property(GenotypeMatrix, "nphase")

def test_mat_format_is_abstract():
    assert_abstract_property(GenotypeMatrix, "mat_format")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mat_asformat_is_abstract():
    assert_abstract_method(GenotypeMatrix, "mat_asformat")

def test_tacount_is_abstract():
    assert_abstract_method(GenotypeMatrix, "tacount")

def test_tafreq_is_abstract():
    assert_abstract_method(GenotypeMatrix, "tafreq")

def test_acount_is_abstract():
    assert_abstract_method(GenotypeMatrix, "acount")

def test_afreq_is_abstract():
    assert_abstract_method(GenotypeMatrix, "afreq")

def test_apoly_is_abstract():
    assert_abstract_method(GenotypeMatrix, "apoly")

def test_maf_is_abstract():
    assert_abstract_method(GenotypeMatrix, "maf")

def test_meh_is_abstract():
    assert_abstract_method(GenotypeMatrix, "meh")

def test_gtcount_is_abstract():
    assert_abstract_method(GenotypeMatrix, "gtcount")

def test_gtfreq_is_abstract():
    assert_abstract_method(GenotypeMatrix, "gtfreq")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GenotypeMatrix_is_concrete():
    assert_concrete_function(check_is_GenotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_GenotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GenotypeMatrix(None, "mat")
