import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import check_is_DenseMolecularCoancestryMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.breed.prot.gt.DenseUnphasedGenotyping import DenseUnphasedGenotyping

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def pgmat_sample(shared_datadir):
    data_path = shared_datadir / "sample.vcf"
    pgmat = DensePhasedGenotypeMatrix.from_vcf(data_path)
    yield pgmat

@pytest.fixture
def gmat_sample(pgmat_sample):
    tmp = DenseUnphasedGenotyping()
    out = tmp.genotype(pgmat_sample)
    yield out

@pytest.fixture
def cmat_pgmat_sample(pgmat_sample):
    out = DenseMolecularCoancestryMatrix.from_gmat(pgmat_sample)
    yield out

@pytest.fixture
def cmat_gmat_sample(gmat_sample):
    out = DenseMolecularCoancestryMatrix.from_gmat(gmat_sample)
    yield out

@pytest.fixture
def cmat_sample(pgmat_sample):
    # data_path = shared_datadir / "Song_2016_phased_chr_1000.vcf"
    out = DenseMolecularCoancestryMatrix.from_gmat(pgmat_sample)
    yield out

@pytest.fixture
def Z_mat_int8():
    Z = numpy.array(
       [[-1,  1,  0,  0, -1, -1,  0, -1,  0, -1, -1,  1, -1,  0,  1, -1],
        [-1,  1, -1, -1,  0,  1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1],
        [ 1,  1,  1,  0,  0,  1,  1, -1, -1,  0,  0,  1, -1, -1, -1, -1],
        [-1,  1,  1, -1,  1,  1,  1,  0,  1,  1,  0,  1,  1, -1, -1,  1],
        [ 1,  1,  1, -1,  0,  0,  0, -1, -1,  1, -1,  0,  1,  0, -1, -1],
        [-1, -1,  1,  0, -1,  0,  1, -1, -1,  0, -1,  1,  0, -1, -1,  1],
        [ 0,  0,  0,  1, -1,  0,  1, -1,  1,  0, -1,  1, -1,  0, -1,  0],
        [-1,  1,  0, -1,  0,  1,  1,  0,  0, -1, -1,  1, -1, -1,  0, -1]], 
        dtype = "int8"
    )
    yield Z

@pytest.fixture
def X_mat_int8(Z_mat_int8):
    yield Z_mat_int8 + 1

@pytest.fixture
def A_mat_float64(Z_mat_int8):
    A = ((1.0/Z_mat_int8.shape[1]) * (Z_mat_int8.dot(Z_mat_int8.T))) + 1.0
    yield A

@pytest.fixture
def cmat_numpy(A_mat_float64):
    yield DenseMolecularCoancestryMatrix(A_mat_float64)

@pytest.fixture
def gmat_numpy(X_mat_int8):
    yield DenseGenotypeMatrix(X_mat_int8, ploidy = 2)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseMolecularCoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DenseMolecularCoancestryMatrix, "__init__")

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_mat_fget(cmat_sample):
    # test matrix properties
    assert numpy.all(cmat_sample >= 0.0)   # completely different individuals have coancestry == 0
    assert numpy.all(cmat_sample <= 2.0)   # identical individuals have coancestry == 2
    # test matrix symmetry
    n = cmat_sample.mat.shape[0]
    for i in range(n):
        for j in range(i,n):
            assert cmat_sample[i,j] == cmat_sample[j,i]

################################################################################
###################### Test concrete method functionality ######################
################################################################################

# from GenotypeMatrix
def test_from_gmat(cmat_numpy, gmat_numpy, A_mat_float64):
    assert numpy.all(cmat_numpy == A_mat_float64)
    assert numpy.all(cmat_numpy == DenseMolecularCoancestryMatrix.from_gmat(gmat_numpy))

def test_from_gmat_pgmat_vs_gmat(cmat_pgmat_sample, cmat_gmat_sample):
    assert numpy.all(cmat_pgmat_sample.mat == cmat_gmat_sample.mat)

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_DenseMolecularCoancestryMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseMolecularCoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseMolecularCoancestryMatrix(cmat_sample):
    with not_raises(TypeError):
        check_is_DenseMolecularCoancestryMatrix(cmat_sample, "cmat_sample")
    with pytest.raises(TypeError):
        check_is_DenseMolecularCoancestryMatrix(None, "cmat_sample")
