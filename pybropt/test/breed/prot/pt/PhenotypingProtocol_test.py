import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.prot.pt import PhenotypingProtocol
from pybropt.breed.prot.pt import is_PhenotypingProtocol
from pybropt.breed.prot.pt import check_is_PhenotypingProtocol
from pybropt.breed.prot.pt import cond_check_is_PhenotypingProtocol


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ptprot():
    yield PhenotypingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(PhenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PhenotypingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_gpmod_is_abstract():
    generic_assert_abstract_property(PhenotypingProtocol, "gpmod")

def test_var_err_is_abstract():
    generic_assert_abstract_property(PhenotypingProtocol, "var_err")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_phenotype_is_abstract():
    generic_assert_abstract_method(PhenotypingProtocol, "phenotype")

def test_set_h2_is_abstract():
    generic_assert_abstract_method(PhenotypingProtocol, "set_h2")

def test_set_H2_is_abstract():
    generic_assert_abstract_method(PhenotypingProtocol, "set_H2")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_PhenotypingProtocol_is_concrete():
    generic_assert_concrete_function(is_PhenotypingProtocol)

def test_check_is_PhenotypingProtocol_is_concrete():
    generic_assert_concrete_function(check_is_PhenotypingProtocol)

def test_cond_check_is_PhenotypingProtocol_is_concrete():
    generic_assert_concrete_function(cond_check_is_PhenotypingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhenotypingProtocol(ptprot):
    assert is_PhenotypingProtocol(ptprot)

def test_check_is_PhenotypingProtocol(ptprot):
    with not_raises(TypeError):
        check_is_PhenotypingProtocol(ptprot, "ptprot")
    with pytest.raises(TypeError):
        check_is_PhenotypingProtocol(None, "ptprot")