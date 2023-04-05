import pytest

from pybrops.test import assert_abstract_methods
from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.breed.prot.pt.PhenotypingProtocol import is_PhenotypingProtocol
from pybrops.breed.prot.pt.PhenotypingProtocol import check_is_PhenotypingProtocol

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
    assert_docstring(PhenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PhenotypingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_gpmod_is_abstract():
    assert_abstract_property(PhenotypingProtocol, "gpmod")

def test_var_err_is_abstract():
    assert_abstract_property(PhenotypingProtocol, "var_err")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_phenotype_is_abstract():
    assert_abstract_method(PhenotypingProtocol, "phenotype")

def test_set_h2_is_abstract():
    assert_abstract_method(PhenotypingProtocol, "set_h2")

def test_set_H2_is_abstract():
    assert_abstract_method(PhenotypingProtocol, "set_H2")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_PhenotypingProtocol_is_concrete():
    assert_concrete_function(is_PhenotypingProtocol)

def test_check_is_PhenotypingProtocol_is_concrete():
    assert_concrete_function(check_is_PhenotypingProtocol)

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
