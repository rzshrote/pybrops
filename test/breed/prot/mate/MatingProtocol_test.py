import numpy
import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function
from .common_fixtures import *
from pybrops.breed.prot.mate.MatingProtocol import MatingProtocol
from pybrops.breed.prot.mate.MatingProtocol import check_is_MatingProtocol

class DummyMatingProtocol(MatingProtocol):
    def __init__(self):
        pass
    def mate(self, pgmat, xconfig, nmating, nprogeny, miscout, **kwargs):
        pass
    @property
    def nparent(self) -> int:
        """nparent."""
        return 0
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set nparent."""
        self._nparent = value
    

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mprot():
    yield DummyMatingProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(MatingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(MatingProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    assert_abstract_method(MatingProtocol, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_MatingProtocol_is_concrete():
    assert_concrete_function(check_is_MatingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MatingProtocol(mprot):
    with not_raises(TypeError):
        check_is_MatingProtocol(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_MatingProtocol(None, "mprot")
