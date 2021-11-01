import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.prot.sel import TwoWayOptimalHaploidValueParentSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat import DensePhasedGenotypeMatrix

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
    yield numpy.int64([1, 1, 1, 1, 1, 2, 2, 2, 2, 2])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_genpos():
    yield numpy.float64([
        0.31, 0.38, 0.66, 0.67, 0.73,
        0.3 , 0.44, 0.58, 0.96, 0.99
    ])

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dpgvmat(mat_int8, mat_chrgrp, mat_phypos, mat_genpos, mat_taxa, mat_taxa_grp):
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

################################################################################
################################ Genomic model #################################
################################################################################
@pytest.fixture
def mu():
    # yield numpy.float64([
    #     [1.4]
    # ])
    yield numpy.float64([
        [1.4],
        [2.5],
        [7.2]
    ])

@pytest.fixture
def beta():
    # yield numpy.float64([
    #     [-0.33],
    #     [-0.69],
    #     [ 1.12],
    #     [-1.44],
    #     [ 0.88],
    #     [ 1.23],
    #     [ 0.19],
    #     [-2.12],
    #     [-0.87],
    #     [ 0.06]
    # ])
    yield numpy.float64([
        [-0.33,  2.08, -2.42],
        [-0.69, -1.87, -1.38],
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
def trait():
    # yield numpy.object_(["protein"])
    yield numpy.object_(["protein", "yield", "quality"])

@pytest.fixture
def model_name():
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(mu, beta, trait, model_name, params):
    yield GenericLinearGenomicModel(
        mu = mu,
        beta = beta,
        trait = trait,
        model_name = model_name,
        params = params
    )

################################################################################
############################ Breeding values model #############################
################################################################################
@pytest.fixture
def bvmat(glgmod, dpgvmat):
    yield glgmod.predict(dpgvmat)

################################################################################
###################### TwoWayOptimalHaploidValueParentSelection ######################
################################################################################
@pytest.fixture
def k_p():
    yield 2

@pytest.fixture
def b_p():
    yield 3

@pytest.fixture
def traitwt_p():
    # yield numpy.float64([1.0])
    yield numpy.float64([1.0, 1.0, 1.0])

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
def twohvps(k_p, traitwt_p, b_p, ncross, nprogeny, rng):
    yield TwoWayOptimalHaploidValueParentSelection(
        k_p = k_p,
        traitwt_p = traitwt_p,
        b_p = b_p,
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################
def test_pselect(twohvps, dpgvmat, bvmat, glgmod, ncross, nprogeny):
    geno = {
        "cand" : dpgvmat,
        "main" : dpgvmat,
        "queue" : [dpgvmat]
    }
    bval = {
        "cand" : bvmat,
        "cand_true" : bvmat,
        "main" : bvmat,
        "main_true" : bvmat
    }
    gmod = {
        "cand" : glgmod,
        "main" : glgmod,
        "true" : glgmod
    }

    out_gmat, out_sel, out_ncross, out_nprogeny, out_misc = twohvps.pselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod
    )

    assert numpy.all(out_sel == [1,3,0,3])
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny
