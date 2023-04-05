import pytest

from pybrops.test import assert_abstract_methods
from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.breed.op.log.Logbook import Logbook
from pybrops.breed.op.log.Logbook import is_Logbook
from pybrops.breed.op.log.Logbook import check_is_Logbook

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def logbook():
    yield Logbook()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(Logbook)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(Logbook, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_data_is_abstract():
    assert_abstract_property(Logbook, "data")

def test_rep_is_abstract():
    assert_abstract_property(Logbook, "rep")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_log_initialize_is_abstract():
    assert_abstract_method(Logbook, "log_initialize")

def test_log_pselect_is_abstract():
    assert_abstract_method(Logbook, "log_pselect")

def test_log_mate_is_abstract():
    assert_abstract_method(Logbook, "log_mate")

def test_log_evaluate_is_abstract():
    assert_abstract_method(Logbook, "log_evaluate")

def test_log_sselect_is_abstract():
    assert_abstract_method(Logbook, "log_sselect")

def test_reset_is_abstract():
    assert_abstract_method(Logbook, "reset")

def test_write_is_abstract():
    assert_abstract_method(Logbook, "write")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_Logbook_is_concrete():
    assert_concrete_function(is_Logbook)

def test_check_is_Logbook_is_concrete():
    assert_concrete_function(check_is_Logbook)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_Logbook(logbook):
    assert is_Logbook(logbook)

def test_check_is_Logbook(logbook):
    with not_raises(TypeError):
        check_is_Logbook(logbook, "logbook")
    with pytest.raises(TypeError):
        check_is_Logbook(None, "logbook")
