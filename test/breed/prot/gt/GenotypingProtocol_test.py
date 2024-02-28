import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol
from pybrops.breed.prot.gt.GenotypingProtocol import check_is_GenotypingProtocol
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def bvprot():
    yield DummyGenotypingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(GenotypingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_genotype_is_abstract():
    assert_method_isabstract(GenotypingProtocol, "genotype")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GenotypingProtocol_is_concrete():
    assert_function_isconcrete(check_is_GenotypingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GenotypingProtocol(bvprot):
    with not_raises(TypeError):
        check_is_GenotypingProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_GenotypingProtocol(None, "bvprot")
