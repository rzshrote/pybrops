import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybrops.model.vmat.DenseGeneticVarianceMatrix import check_is_DenseGeneticVarianceMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.array([
        [[1.5, 2.3, 8.3],
         [2.3, 0.4, 3.6],
         [8.3, 3.6, 0.9]],
        [[0.7, 1.1, 4.8],
         [1.1, 0.2, 3.2],
         [4.8, 3.2, 0.4]]
    ])
    a = a.transpose(1,2,0)
    yield a

@pytest.fixture
def mat_taxa():
    a = numpy.object_(["A", "B", "C"])
    yield a

@pytest.fixture
def mat_taxa_grp():
    a = numpy.int64([0,1,1])
    yield a

@pytest.fixture
def mat(mat_float64, mat_taxa, mat_taxa_grp):
    yield DummyDenseGeneticVarianceMatrix(
        mat = mat_float64,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseGeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseGeneticVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseGeneticVarianceMatrix_is_concrete():
    assert_concrete_function(check_is_DenseGeneticVarianceMatrix)

def test_check_is_DenseGeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseGeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseGeneticVarianceMatrix(None, "mat")