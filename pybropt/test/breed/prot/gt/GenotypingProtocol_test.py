import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.prot.gt import GenotypingProtocol
from pybropt.breed.prot.gt import is_GenotypingProtocol
from pybropt.breed.prot.gt import check_is_GenotypingProtocol
from pybropt.breed.prot.gt import cond_check_is_GenotypingProtocol


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def bvprot():
    yield GenotypingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GenotypingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_genotype_is_abstract():
    generic_assert_abstract_method(GenotypingProtocol, "genotype")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_GenotypingProtocol_is_concrete():
    generic_assert_concrete_function(is_GenotypingProtocol)

def test_check_is_GenotypingProtocol_is_concrete():
    generic_assert_concrete_function(check_is_GenotypingProtocol)

def test_cond_check_is_GenotypingProtocol_is_concrete():
    generic_assert_concrete_function(cond_check_is_GenotypingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GenotypingProtocol(bvprot):
    assert is_GenotypingProtocol(bvprot)

def test_check_is_GenotypingProtocol(bvprot):
    with not_raises(TypeError):
        check_is_GenotypingProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_GenotypingProtocol(None, "bvprot")