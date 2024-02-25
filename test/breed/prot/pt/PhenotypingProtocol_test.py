import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.breed.prot.pt.PhenotypingProtocol import check_is_PhenotypingProtocol
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ptprot():
    yield DummyPhenotypingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(PhenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(PhenotypingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_gpmod_is_abstract():
    assert_property_isabstract(PhenotypingProtocol, "gpmod")

def test_var_err_is_abstract():
    assert_property_isabstract(PhenotypingProtocol, "var_err")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_phenotype_is_abstract():
    assert_method_isabstract(PhenotypingProtocol, "phenotype")

def test_set_h2_is_abstract():
    assert_method_isabstract(PhenotypingProtocol, "set_h2")

def test_set_H2_is_abstract():
    assert_method_isabstract(PhenotypingProtocol, "set_H2")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_PhenotypingProtocol_is_concrete():
    assert_function_isconcrete(check_is_PhenotypingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhenotypingProtocol(ptprot):
    with not_raises(TypeError):
        check_is_PhenotypingProtocol(ptprot, "ptprot")
    with pytest.raises(TypeError):
        check_is_PhenotypingProtocol(None, "ptprot")
