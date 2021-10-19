import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.prot.sel import OptimalContributionSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.cmat import DenseMolecularCoancestryMatrix

################################################################################
################################## Genotypes ###################################
################################################################################
@pytest.fixture
def mat_int8():
    yield numpy.int8([
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
def dpgvmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DensePhasedGenotypeMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

################################################################################
################################ Genomic model #################################
################################################################################
@pytest.fixture
def beta():
    yield numpy.float64([[1.4]])
    # yield numpy.float64([[1.4, 2.5, 7.2]])

@pytest.fixture
def u():
    yield numpy.float64([
        [-0.33],
        [-0.69],
        [ 1.12],
        [-1.44],
        [ 0.88],
        [ 1.23],
        [ 0.19],
        [-2.12],
        [-0.87],
        [ 0.06]
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
    yield numpy.object_(["protein"])
    # yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(beta, u, trait, model_name, params):
    yield GenericLinearGenomicModel(
        beta = beta,
        u = u,
        trait = trait,
        model_name = model_name,
        params = params
    )

################################################################################
############################ Breeding values model #############################
################################################################################
@pytest.fixture
def bvmat(glgmod, dpgvmat):
    yield glgmod.gebv(dpgvmat)

################################################################################
###################### OptimalContributionSelection ######################
################################################################################
@pytest.fixture
def nparent():
    yield 2

@pytest.fixture
def traitwt_p():
    yield numpy.float64([1.0])
    # yield numpy.float64([1.0, 1.0, 1.0])

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def rng():
    yield Generator(PCG64(192837465))

@pytest.fixture
def inbfn():
    def fn(t_cur, t_max):
        return 0.75
    yield fn

@pytest.fixture
def cmatcls():
    yield DenseMolecularCoancestryMatrix

@pytest.fixture
def ocs(nparent, traitwt_p, inbfn, ncross, nprogeny, cmatcls, rng):
    yield OptimalContributionSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        traitwt_p = traitwt_p,
        inbfn = inbfn,
        cmatcls = cmatcls,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################
