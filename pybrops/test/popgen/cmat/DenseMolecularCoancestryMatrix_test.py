import pytest
import numpy

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import is_DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import check_is_DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import cond_check_is_DenseMolecularCoancestryMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def cmat(shared_datadir):
    # data_path = shared_datadir / "Song_2016_phased_chr_1000.vcf"
    data_path = shared_datadir / "sample.vcf"
    gmat = DensePhasedGenotypeMatrix.from_vcf(data_path)
    out = DenseMolecularCoancestryMatrix.from_gmat(gmat)
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseMolecularCoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseMolecularCoancestryMatrix, "__init__")

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_mat_fget(cmat):
    # test matrix properties
    assert numpy.all(cmat >= 0.0)   # completely different individuals
    assert numpy.all(cmat <= 1.0)   # identical individuals
    n = cmat.mat.shape[0]
    for i in range(n):
        for j in range(i,n):
            assert cmat[i,j] == cmat[j,i]

################################################################################
###################### Test concrete method functionality ######################
################################################################################

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DenseMolecularCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseMolecularCoancestryMatrix)

def test_check_is_DenseMolecularCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseMolecularCoancestryMatrix)

def test_cond_check_is_DenseMolecularCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseMolecularCoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseMolecularCoancestryMatrix(cmat):
    assert is_DenseMolecularCoancestryMatrix(cmat)

def test_check_is_DenseMolecularCoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_DenseMolecularCoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_DenseMolecularCoancestryMatrix(None, "cmat")
