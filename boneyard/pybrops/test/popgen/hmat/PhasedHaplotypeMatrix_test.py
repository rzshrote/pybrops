import inspect
import pytest

from pybrops.test import generic_test_abstract_methods

from pybrops.popgen.hmat import PhasedHaplotypeMatrix
from pybrops.popgen.hmat import is_PhasedHaplotypeMatrix
from pybrops.popgen.hmat import check_is_PhasedHaplotypeMatrix
from pybrops.popgen.hmat import cond_check_is_PhasedHaplotypeMatrix

@pytest.fixture
def phmat():
    yield PhasedHaplotypeMatrix()

@pytest.fixture
def vmethods(phmat):
    yield [m for m in dir(phmat) if m.startswith('__') is False]

################################################################################
#################################### Tests #####################################
################################################################################
def test_abstract_methods(phmat, vmethods):
    generic_test_abstract_methods(phmat, vmethods)

def test_is_PhasedHaplotypeMatrix(phmat):
    assert is_PhasedHaplotypeMatrix(phmat)

def test_check_is_PhasedHaplotypeMatrix(phmat):
    check_is_PhasedHaplotypeMatrix(phmat, "phmat")
    with pytest.raises(TypeError):
        check_is_PhasedHaplotypeMatrix(None, "None")

def test_cond_check_is_PhasedHaplotypeMatrix(phmat):
    cond_check_is_PhasedHaplotypeMatrix(phmat, "phmat")
    cond_check_is_PhasedHaplotypeMatrix(None, "None")
    with pytest.raises(TypeError):
        cond_check_is_PhasedHaplotypeMatrix("incorrect input type", "str")
