import numpy
import pytest
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.fcty.DenseAdditiveGenicVarianceMatrixFactory import DenseAdditiveGenicVarianceMatrixFactory
from pybrops.model.vmat.fcty.DenseAdditiveGenicVarianceMatrixFactory import check_is_DenseAdditiveGenicVarianceMatrixFactory

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def mat_beta():
    yield numpy.float64([
        [1.4, 2.5, 7.2]
    ])

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a():
    yield numpy.float64([
        [-0.33,  2.08, -2.42],
        [-0.69, -1.87,  1.38],
        [ 1.12,  1.38, -5.65],
        [-1.44,  0.20,  4.22],
        [ 0.88, -0.81,  1.55],
        [ 1.23,  0.25,  5.13],
        [ 0.19,  4.35,  0.15],
        [-2.12,  0.73, -0.38],
        [-0.87,  1.25,  2.38],
        [ 0.06, -2.52,  2.48]
    ])

@pytest.fixture
def mat_trait():
    yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_dalgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def dalgmod(mat_beta, mat_u_misc, mat_u_a, mat_trait, model_name, params):
    yield DenseAdditiveLinearGenomicModel(
        beta = mat_beta,
        u_misc = mat_u_misc,
        u_a = mat_u_a,
        trait = mat_trait,
        model_name = model_name,
        params = params
    )

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[1, 0, 1, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 0, 0, 1, 1],
        [1, 0, 1, 0, 0, 1, 0, 0, 0, 1]],
       [[0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

@pytest.fixture
def mat_genpos():
    yield numpy.array([0.5, 0.7, 0.1, 0.4, 0.2, 0.9, 0.4, 0.7, 0.2, 0.6], dtype = float)

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dpgmat(mat_int8, mat_chrgrp, mat_phypos, mat_genpos, mat_taxa, mat_taxa_grp):
    out = DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        vrnt_genpos = mat_genpos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )
    out.group_vrnt()
    out.group_taxa()
    yield out

############################################################
##################### Variance Matrix ######################
############################################################

############################################################
######################### Factory ##########################
############################################################
@pytest.fixture
def fcty():
    yield DenseAdditiveGenicVarianceMatrixFactory()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseAdditiveGenicVarianceMatrixFactory)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseAdditiveGenicVarianceMatrixFactory, "__init__")

def test_from_gmod_is_concrete():
    assert_concrete_method(DenseAdditiveGenicVarianceMatrixFactory, "from_gmod")

def test_from_algmod_is_concrete():
    assert_concrete_method(DenseAdditiveGenicVarianceMatrixFactory, "from_algmod")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################## Test object methods #############################
################################################################################
def test_from_gmod(fcty, dalgmod, dpgmat, mat_taxa, mat_taxa_grp):
    vmat = fcty.from_gmod(
        gmod = dalgmod,
        pgmat = dpgmat,
        nprogeny = 40,
    )
    assert isinstance(vmat, DenseAdditiveGenicVarianceMatrix)
    assert numpy.all(vmat.taxa == mat_taxa)
    assert numpy.all(vmat.taxa_grp == mat_taxa_grp)

def test_from_algmod(fcty, dalgmod, dpgmat, mat_taxa, mat_taxa_grp):
    vmat = fcty.from_algmod(
        algmod = dalgmod,
        pgmat = dpgmat,
        nprogeny = 40,
    )
    assert isinstance(vmat, DenseAdditiveGenicVarianceMatrix)
    assert numpy.all(vmat.taxa == mat_taxa)
    assert numpy.all(vmat.taxa_grp == mat_taxa_grp)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseAdditiveGenicVarianceMatrixFactory_is_concrete():
    assert_concrete_function(check_is_DenseAdditiveGenicVarianceMatrixFactory)

def test_check_is_DenseAdditiveGenicVarianceMatrixFactory(fcty):
    with not_raises(TypeError):
        check_is_DenseAdditiveGenicVarianceMatrixFactory(fcty, "fcty")
    with pytest.raises(TypeError):
        check_is_DenseAdditiveGenicVarianceMatrixFactory(None, "fcty")
