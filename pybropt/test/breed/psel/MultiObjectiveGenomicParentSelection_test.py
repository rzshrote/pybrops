import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.psel import MultiObjectiveGenomicParentSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.algo.opt import SteepestAscentSetHillClimber
from pybropt.algo.opt import StochasticAscentSetHillClimber

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
def mat_int8_big():
    yield numpy.int8([[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]],

       [[1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]]])

@pytest.fixture
def mat_chrgrp():
    yield numpy.int64([1, 1, 2, 2, 3, 3, 4, 4, 5, 5])

@pytest.fixture
def mat_chrgrp_big():
    yield numpy.int64([1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2])

@pytest.fixture
def mat_phypos():
    yield numpy.int64([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

@pytest.fixture
def mat_phypos_big():
    yield numpy.arange(16)

@pytest.fixture
def mat_taxa():
    yield numpy.object_(["Line1", "Line2", "Line3", "Line4", "Line5"])

@pytest.fixture
def mat_taxa_big():
    yield numpy.object_(["Line"+str(i).zfill(2) for i in range(30)])

@pytest.fixture
def mat_taxa_grp():
    yield numpy.int64([1, 1, 2, 2, 2])

@pytest.fixture
def dpgvmat(mat_int8, mat_chrgrp, mat_phypos, mat_taxa, mat_taxa_grp):
    yield DensePhasedGenotypeVariantMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

@pytest.fixture
def dpgvmat_big(mat_int8_big, mat_chrgrp_big, mat_phypos_big, mat_taxa_big):
    yield DensePhasedGenotypeVariantMatrix(
        mat = mat_int8_big,
        vrnt_chrgrp = mat_chrgrp_big,
        vrnt_phypos = mat_phypos_big,
        taxa = mat_taxa_big
    )

################################################################################
################################ Genomic model #################################
################################################################################
@pytest.fixture
def mu():
    yield numpy.float64([
        [1.4],
        [2.5],
        [7.2]
    ])

@pytest.fixture
def beta():
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
def beta_big():
    yield numpy.float64([
       [-0.87, -0.16, -0.04],
       [-0.03,  0.05, -0.15],
       [ 0.36, -0.15,  0.54],
       [ 2.35,  1.12,  0.33],
       [-0.93, -0.59, -1.1 ],
       [-0.63,  0.61,  0.15],
       [-0.8 ,  0.95, -0.56],
       [-1.03,  1.55,  0.12],
       [-1.21,  0.79,  1.42],
       [-0.09, -1.88, -1.83],
       [-0.03,  1.97, -1.98],
       [-0.04,  0.2 ,  1.43],
       [ 0.45,  1.26, -2.21],
       [-0.31, -0.62,  1.09],
       [ 0.9 , -1.37,  0.91],
       [-0.87,  1.2 , -1.68]
   ])

@pytest.fixture
def trait():
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

@pytest.fixture
def glgmod_big(mu, beta_big, trait, model_name, params):
    yield GenericLinearGenomicModel(
        mu = mu,
        beta = beta_big,
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

@pytest.fixture
def bvmat_big(glgmod_big, dpgvmat_big):
    yield glgmod_big.predict(dpgvmat_big)

################################################################################
###################### MultiObjectiveGenomicParentSelection ######################
################################################################################
@pytest.fixture
def k_p():
    yield 2

@pytest.fixture
def traitobjwt_p():
    yield numpy.float64([[1.0, 1.0, 1.0]])

@pytest.fixture
def traitsum_p():
    yield True

@pytest.fixture
def objsum_p():
    yield True

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
def algorithm_p(k_p, dpgvmat, rng):
    yield StochasticAscentSetHillClimber(
        k = k_p,
        setspace = numpy.arange(dpgvmat.ntaxa),
        rng = rng,
        objwt = 1.0,
    )
    # yield SteepestAscentSetHillClimber(
    #     k = k_p,
    #     setspace = numpy.arange(dpgvmat.ntaxa),
    #     rng = rng,
    #     objwt = 1.0,
    # )

@pytest.fixture
def mogps(k_p, traitobjwt_p, traitsum_p, objsum_p, algorithm_p, ncross, nprogeny, rng):
    yield MultiObjectiveGenomicParentSelection(
        k_p = k_p,
        traitobjwt_p = traitobjwt_p,
        traitsum_p = traitsum_p,
        objsum_p = objsum_p,
        algorithm_p = algorithm_p,
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################
def test_pselect(mogps, dpgvmat, bvmat, glgmod, ncross, nprogeny):
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

    out_gmat, out_sel, out_ncross, out_nprogeny, out_misc = mogps.pselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod
    )

    assert numpy.all(out_sel == [1,2]) or numpy.all(out_sel == [2,1])
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny

def test_ppareto(mogps, dpgvmat_big, bvmat_big, glgmod_big, ncross, nprogeny):
    geno = {
        "cand" : dpgvmat_big,
        "main" : dpgvmat_big,
        "queue" : [dpgvmat_big]
    }
    bval = {
        "cand" : bvmat_big,
        "cand_true" : bvmat_big,
        "main" : bvmat_big,
        "main_true" : bvmat_big
    }
    gmod = {
        "cand" : glgmod_big,
        "main" : glgmod_big,
        "true" : glgmod_big
    }

    frontier, pop, logbook = mogps.ppareto(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod,
        k = 6
    )

    assert isinstance(frontier, numpy.ndarray)

    from matplotlib import pyplot

    pyplot.scatter(frontier[:,0], frontier[:,1], c="b")
    pyplot.axis("tight")
    pyplot.savefig("frontier.png")

    # for i in range(len(pop)):
    #     print(pop[i])
    # raise RuntimeError
