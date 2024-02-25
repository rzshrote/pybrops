import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.op.log.Logbook import Logbook
from pybrops.breed.op.log.Logbook import check_is_Logbook
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def logbook():
    yield DummyLogbook()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(Logbook)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(Logbook, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_data_is_abstract():
    assert_property_isabstract(Logbook, "data")

def test_rep_is_abstract():
    assert_property_isabstract(Logbook, "rep")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_log_initialize_is_abstract():
    assert_method_isabstract(Logbook, "log_initialize")

def test_log_pselect_is_abstract():
    assert_method_isabstract(Logbook, "log_pselect")

def test_log_mate_is_abstract():
    assert_method_isabstract(Logbook, "log_mate")

def test_log_evaluate_is_abstract():
    assert_method_isabstract(Logbook, "log_evaluate")

def test_log_sselect_is_abstract():
    assert_method_isabstract(Logbook, "log_sselect")

def test_reset_is_abstract():
    assert_method_isabstract(Logbook, "reset")

def test_write_is_abstract():
    assert_method_isabstract(Logbook, "write")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_Logbook_is_concrete():
    assert_function_isconcrete(check_is_Logbook)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_Logbook(logbook):
    with not_raises(TypeError):
        check_is_Logbook(logbook, "logbook")
    with pytest.raises(TypeError):
        check_is_Logbook(None, "logbook")
