import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol
from pybrops.breed.prot.gt.GenotypingProtocol import check_is_GenotypingProtocol

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
    assert_docstring(GenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GenotypingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_genotype_is_abstract():
    assert_abstract_method(GenotypingProtocol, "genotype")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GenotypingProtocol_is_concrete():
    assert_concrete_function(check_is_GenotypingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenotypingProtocol(bvprot):
    with not_raises(TypeError):
        check_is_GenotypingProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_GenotypingProtocol(None, "bvprot")
