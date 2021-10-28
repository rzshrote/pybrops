import numpy
import pytest
from numpy.random import Generator
from numpy.random import PCG64
from matplotlib import pyplot

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.prot.pt import TruePhenotyping
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.ptdf import is_PhenotypeDataFrame

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    a = numpy.int8([
       [[1, 0, 0, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 1]],
       [[0, 0, 1, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
    ])
    yield a

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

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
def dpgmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def beta():
    yield numpy.float64([[1.4, 2.5]])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def u():
    yield numpy.float64([
        [-0.33,  2.08],
        [-0.69, -1.87],
        [ 1.12,  1.38],
        [-1.44,  0.20],
        [ 0.88, -0.81],
        [ 1.23,  0.25],
        [ 0.19,  4.35],
        [-2.12,  0.73],
        [-0.87,  1.25],
        [ 0.06, -2.52]
    ])
    # yield numpy.float64([
    #     [-0.33,  2.08, -2.42],
    #     [-0.69, -1.87, -1.38],
    #     [ 1.12,  1.38, -5.65],
    #     [-1.44,  0.20,  4.22],
    #     [ 0.88, -0.81,  1.55],
    #     [ 1.23,  0.25,  5.13],
    #     [ 0.19,  4.35,  0.15],
    #     [-2.12,  0.73, -0.38],
    #     [-0.87,  1.25,  2.38],
    #     [ 0.06, -2.52,  2.48]
    # ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield"])
    # yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_gpmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def gpmod(beta, u, trait, model_name, params):
    yield GenericLinearGenomicModel(
        beta = beta,
        u = u,
        trait = trait,
        model_name = model_name,
        params = params
    )

############################################################
################## Breeding values model ###################
############################################################
@pytest.fixture
def bvmat(gpmod, dpgmat):
    yield gpmod.gebv(dpgmat)

############################################################
##################### TruePhenotyping ######################
############################################################
@pytest.fixture
def ptprot():
    yield TruePhenotyping()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(TruePhenotyping)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_phenotype_is_concrete():
    generic_assert_concrete_method(TruePhenotyping, "phenotype")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_phenotype(ptprot, dpgmat, gpmod):
    df = ptprot.phenotype(dpgmat, gpmod)
    assert is_PhenotypeDataFrame(df)

    expected_ncol = gpmod.ntrait
    if dpgmat.taxa is not None:
        expected_ncol += 1
    if dpgmat.taxa_grp is not None:
        expected_ncol += 1
    assert df.ncol == expected_ncol

    expected_nrow = dpgmat.ntaxa
    assert df.nrow == expected_nrow
