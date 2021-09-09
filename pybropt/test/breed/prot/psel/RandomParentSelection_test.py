import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.psel import RandomParentSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix

################################################################################
################################## Genotypes ###################################
################################################################################
@pytest.fixture
def mat_int8_big():
    yield numpy.int8([
       [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
        [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]]
    ])

@pytest.fixture
def mat_chrgrp_big():
    yield numpy.int64([1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2])

@pytest.fixture
def mat_phypos_big():
    yield numpy.arange(16)

@pytest.fixture
def mat_taxa_big():
    yield numpy.object_(["Line"+str(i).zfill(2) for i in range(30)])

@pytest.fixture
def mat_taxa_grp_big():
    yield numpy.int64([
        1,1,1,1,1,1,
        2,2,2,2,2,2,
        3,3,3,3,3,3,
        4,4,4,4,4,4,
        5,5,5,5,5,5
    ])

@pytest.fixture
def dpgvmat_big(mat_int8_big, mat_chrgrp_big, mat_phypos_big, mat_taxa_big, mat_taxa_grp_big):
    yield DensePhasedGenotypeVariantMatrix(
        mat = mat_int8_big,
        vrnt_chrgrp = mat_chrgrp_big,
        vrnt_phypos = mat_phypos_big,
        taxa = mat_taxa_big,
        taxa_grp = mat_taxa_grp_big
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
def bvmat_big(glgmod_big, dpgvmat_big):
    yield glgmod_big.predict(dpgvmat_big)

################################################################################
##################### MultiObjectiveGenomicParentSelection #####################
################################################################################
@pytest.fixture
def nparent():
    yield 10

@pytest.fixture
def ncross():
    yield 1

@pytest.fixture
def nprogeny():
    yield 10

@pytest.fixture
def rng():
    yield Generator(PCG64(239987540))

@pytest.fixture
def rps(nparent, ncross, nprogeny, rng):
    yield RandomParentSelection(
        nparent = nparent,
        ncross = ncross,
        nprogeny = nprogeny,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################
# test constructor
def test_init(rps):
    assert True

# test instance data
def test_nparent(rps, nparent):
    assert rps.nparent == nparent

def test_ncross(rps, ncross):
    assert rps.ncross == ncross

def test_nprogeny(rps, nprogeny):
    assert rps.nprogeny == nprogeny

# test selection criteria
def test_pselect(rps, dpgvmat_big, bvmat_big, glgmod_big, nparent, ncross, nprogeny):
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

    out_gmat, out_sel, out_ncross, out_nprogeny, out_misc = rps.pselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod
    )

    assert len(out_sel) == nparent
    assert len(numpy.unique(out_sel)) == nparent
    assert out_ncross == ncross
    assert out_nprogeny == nprogeny

def test_pobjfn(rps):
    with pytest.raises(RuntimeError):
        rps.pobjfn(None, None, None, None, None, None, None)

def test_pobjfn_vec(rps):
    with pytest.raises(RuntimeError):
        rps.pobjfn_vec(None, None, None, None, None)

def test_ppareto(rps):
    with pytest.raises(RuntimeError):
        rps.ppareto(None, None, None, None, None)
