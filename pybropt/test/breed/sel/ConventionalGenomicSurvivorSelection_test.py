import numpy
import pytest

from numpy.random import Generator
from numpy.random import PCG64

from pybropt.breed.sel import ConventionalGenomicSurvivorSelection
from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix

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
    yield numpy.string_(["Line1", "Line2", "Line3", "Line4", "Line5"])

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
def trait():
    yield numpy.string_(["protein", "yield", "quality"])

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
###################### ConventionalGenomicSurvivorSelection ######################
################################################################################
@pytest.fixture
def k_s():
    yield 2

@pytest.fixture
def traitwt_s():
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
def cgps(k_s, traitwt_s, ncross, nprogeny, rng):
    yield ConventionalGenomicSurvivorSelection(
        k_s = k_s,
        traitwt_s = traitwt_s,
        rng = rng
    )

################################################################################
#################################### Tests #####################################
################################################################################
def test_sselect(cgps, dpgvmat, bvmat, glgmod, ncross, nprogeny):
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

    geno_dict, bval_dict, gmod_dict, misc = cgps.sselect(
        t_cur = 0,
        t_max = 20,
        geno = geno,
        bval = bval,
        gmod = gmod
    )

    assert isinstance(geno_dict, dict)
    assert "cand" in geno_dict
    assert "main" in geno_dict
    assert "queue" in geno_dict

    assert isinstance(bval_dict, dict)
    assert "cand" in bval_dict
    assert "cand_true" in bval_dict
    assert "main" in bval_dict
    assert "main_true" in bval_dict

    assert isinstance(gmod_dict, dict)
    assert "cand" in gmod_dict
    assert "main" in gmod_dict
    assert "true" in gmod_dict
