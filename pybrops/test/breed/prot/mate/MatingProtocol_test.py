import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.breed.prot.mate.MatingProtocol import is_MatingProtocol
from pybrops.breed.prot.mate.MatingProtocol import check_is_MatingProtocol
from pybrops.breed.prot.mate.MatingProtocol import cond_check_is_MatingProtocol


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mprot():
    yield MatingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(MatingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(MatingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    generic_assert_abstract_method(MatingProtocol, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_MatingProtocol_is_concrete():
    generic_assert_concrete_function(is_MatingProtocol)

def test_check_is_MatingProtocol_is_concrete():
    generic_assert_concrete_function(check_is_MatingProtocol)

def test_cond_check_is_MatingProtocol_is_concrete():
    generic_assert_concrete_function(cond_check_is_MatingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_MatingProtocol(mprot):
    assert is_MatingProtocol(mprot)

def test_check_is_MatingProtocol(mprot):
    with not_raises(TypeError):
        check_is_MatingProtocol(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_MatingProtocol(None, "mprot")
