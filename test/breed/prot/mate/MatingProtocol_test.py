import numpy
import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
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
    assert_class_documentation(MatingProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    assert_method_isabstract(MatingProtocol, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_MatingProtocol_is_concrete():
    assert_function_isconcrete(check_is_MatingProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MatingProtocol(mprot):
    with not_raises(TypeError):
        check_is_MatingProtocol(mprot, "mprot")
    with pytest.raises(TypeError):
        check_is_MatingProtocol(None, "mprot")
