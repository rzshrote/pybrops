import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.popgen.hmat import HaplotypeMatrix
from pybropt.popgen.hmat import is_HaplotypeMatrix
from pybropt.popgen.hmat import check_is_HaplotypeMatrix
from pybropt.popgen.hmat import cond_check_is_HaplotypeMatrix

@pytest.fixture
def hmat():
    yield HaplotypeMatrix()

@pytest.fixture
def vmethods(hmat):
    yield [m for m in dir(hmat) if m.startswith('__') is False]

################################################################################
#################################### Tests #####################################
################################################################################
def test_abstract_methods(hmat, vmethods):
    generic_test_abstract_methods(hmat, vmethods)

def test_is_HaplotypeMatrix(hmat):
    assert is_HaplotypeMatrix(hmat)

def test_check_is_HaplotypeMatrix(hmat):
    check_is_HaplotypeMatrix(hmat, "hmat")
    with pytest.raises(TypeError):
        check_is_HaplotypeMatrix(None, "None")

def test_cond_check_is_HaplotypeMatrix(hmat):
    cond_check_is_HaplotypeMatrix(hmat, "hmat")
    cond_check_is_HaplotypeMatrix(None, "None")
    with pytest.raises(TypeError):
        cond_check_is_HaplotypeMatrix("incorrect input type", "str")
