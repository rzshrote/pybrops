import pytest
import numpy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import DenseTwoWayDHAdditiveGeneticVarianceMatrix
# from pybropt.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import is_DenseTwoWayDHAdditiveGeneticVarianceMatrix
# from pybropt.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import check_is_DenseTwoWayDHAdditiveGeneticVarianceMatrix
# from pybropt.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import cond_check_is_DenseTwoWayDHAdditiveGeneticVarianceMatrix

from pybropt.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybropt.popgen.gmap.HaldaneMapFunction import HaldaneMapFunction

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
############### Minimum object construction ################
############################################################
@pytest.fixture
def mat_float64_min():
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
def mat_taxa_min():
    a = numpy.object_(["A", "B", "C"])
    yield a

@pytest.fixture
def mat_taxa_grp_min():
    a = numpy.int64([0,1,1])
    yield a

@pytest.fixture
def mat(mat_float64_min, mat_taxa_min, mat_taxa_grp_min):
    yield DenseTwoWayDHAdditiveGeneticVarianceMatrix(
        mat = mat_float64_min,
        taxa = mat_taxa_min,
        taxa_grp = mat_taxa_grp_min
    )

############################################################
######################## Genotypes #########################
############################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
       [[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]],

       [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2
    ])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([
         1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20
    ])

@pytest.fixture
def mat_genpos():
    yield numpy.float64([
        0.13, 0.32, 0.53, 0.54, 0.55, 0.61, 0.63, 0.7 , 0.75, 0.96,
        0.14, 0.16, 0.26, 0.31, 0.31, 0.68, 0.7 , 0.74, 0.75, 0.91
    ])

@pytest.fixture
def mat_taxa():
    yield numpy.object_([
        'Line01', 'Line02', 'Line03', 'Line04', 'Line05',
        'Line06', 'Line07', 'Line08', 'Line09', 'Line10',
        'Line11', 'Line12', 'Line13', 'Line14', 'Line15',
        'Line16', 'Line17', 'Line18', 'Line19', 'Line20'
    ])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([
        1, 1, 1, 1, 1,
        2, 2, 2, 2, 2,
        3, 3, 3, 3, 3,
        4, 4, 4, 4, 4
    ])

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
    out.group()
    yield out

############################################################
###################### Genomic model #######################
############################################################
@pytest.fixture
def mat_beta():
    yield numpy.float64([
        [25.6, 13.4]
    ])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def mat_u_misc():
    yield None

@pytest.fixture
def mat_u_a():
    yield numpy.float64([
        [ 1.25, -0.68],
        [-0.02, -1.09],
        [ 0.21, -0.5 ],
        [-2.84,  0.64],
        [-1.37, -0.81],
        [-2.06,  2.22],
        [ 1.52, -0.21],
        [-0.23, -1.78],
        [ 1.04, -0.55],
        [-0.77, -1.4 ],
        [-0.44,  0.89],
        [ 0.12, -0.87],
        [-0.44, -0.55],
        [ 1.36,  0.73],
        [ 1.04,  1.22],
        [-0.05,  0.82],
        [ 0.93,  0.73],
        [-0.89,  1.21],
        [ 0.05, -1.19],
        [-1.27, -2.  ]
    ])

@pytest.fixture
def trait():
    yield numpy.object_(["protein", "yield"])

@pytest.fixture
def model_name():
    yield "test_dalgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def dalgmod(mat_beta, mat_u_misc, mat_u_a, trait, model_name, params):
    yield DenseAdditiveLinearGenomicModel(
        beta = mat_beta,
        u_misc = mat_u_misc,
        u_a = mat_u_a,
        trait = trait,
        model_name = model_name,
        params = params
    )

############################################################
############### Genomic variance estimation ################
############################################################
@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def gmapfn():
    yield HaldaneMapFunction()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseTwoWayDHAdditiveGeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseTwoWayDHAdditiveGeneticVarianceMatrix, "__init__")

# def test_from_gmod_is_concrete():
#     generic_assert_concrete_method(DenseTwoWayDHAdditiveGeneticVarianceMatrix, "from_gmod")

def test_from_algmod_is_concrete():
    generic_assert_concrete_method(DenseTwoWayDHAdditiveGeneticVarianceMatrix, "from_algmod")

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_from_algmod(dalgmod, dpgmat, ncross, nprogeny, gmapfn):
    s0 = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
        algmod = dalgmod,
        pgmat = dpgmat,
        ncross = ncross,
        nprogeny = nprogeny,
        s = 0,
        gmapfn = gmapfn,
        mem = None
    )
    s1 = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
        algmod = dalgmod,
        pgmat = dpgmat,
        ncross = ncross,
        nprogeny = nprogeny,
        s = 1,
        gmapfn = gmapfn,
        mem = None
    )

    assert not numpy.all(s0 == s1)
    for i in range(s0.mat.shape[2]):
        s0s = s0.mat[:,:,i]
        s1s = s1.mat[:,:,i]
        assert numpy.all(s0s == s0s.T)
        assert numpy.all(s1s == s1s.T)
