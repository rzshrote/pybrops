import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import check_is_BreedingValueMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyBreedingValueMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(BreedingValueMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_location_is_abstract():
    assert_property_isabstract(BreedingValueMatrix, "location")

def test_scale_is_abstract():
    assert_property_isabstract(BreedingValueMatrix, "scale")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_targmax_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "targmax")

def test_targmin_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "targmin")

def test_tmax_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "tmax")

def test_tmean_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "tmean")

def test_tmin_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "tmin")

def test_trange_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "trange")

def test_tstd_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "tstd")

def test_tvar_is_abstract():
    assert_method_isabstract(BreedingValueMatrix, "tvar")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_BreedingValueMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingValueMatrix(mat):
    with not_raises(TypeError):
        check_is_BreedingValueMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_BreedingValueMatrix(None, "mat")
