import numpy
import pytest
import os.path

from pybropt.model.gmod import GenericLinearGenomicModel
from pybropt.model.gmod import is_GenericLinearGenomicModel
from pybropt.model.gmod import check_is_GenericLinearGenomicModel
from pybropt.model.gmod import cond_check_is_GenericLinearGenomicModel

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.popgen.bvmat import is_BreedingValueMatrix

################################################################################
############################## Fixtures for tests ##############################
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
    yield DensePhasedGenotypeVariantMatrix(
        mat = mat_int8,
        vrnt_chrgrp = mat_chrgrp,
        vrnt_phypos = mat_phypos,
        taxa = mat_taxa,
        taxa_grp = mat_taxa_grp
    )

################################################################################
#################################### Tests #####################################
################################################################################

########################################
############ Utility tests #############
########################################
def test_is_GenericLinearGenomicModel(glgmod):
    assert is_GenericLinearGenomicModel(glgmod)

def test_check_is_GenericLinearGenomicModel():
    with pytest.raises(TypeError):
        check_is_GenericLinearGenomicModel(None, "None")

def test_cond_check_is_GenericLinearGenomicModel():
    with pytest.raises(TypeError):
        cond_check_is_GenericLinearGenomicModel(0, "0")

########################################
############ Property tests ############
########################################
def test_mu_fget(glgmod, mu):
    assert numpy.all(glgmod.mu == mu)

def test_beta_fget(glgmod, beta):
    assert numpy.all(glgmod.beta == beta)

def test_trait_fget(glgmod, trait):
    assert numpy.all(glgmod.trait == trait)

def test_model_name_fget(glgmod, model_name):
    assert glgmod.model_name == model_name

def test_params_fget(glgmod, params):
    assert glgmod.params == params

################################################################################
################################# Method tests #################################
################################################################################

########################################
########### Prediction tests ###########
########################################
def test_fit(glgmod):
    with pytest.raises(RuntimeError):
        glgmod.fit(None, None)

def test_predict_type(glgmod, dpgvmat):
    assert is_BreedingValueMatrix(glgmod.predict(dpgvmat))

def test_predict_value(glgmod, dpgvmat, mu, beta):
    ghat = mu.T + (dpgvmat.tacount() @ beta)
    assert numpy.all(glgmod.predict(dpgvmat) == ghat)

def test_score_value(glgmod, dpgvmat):
    bvmat = glgmod.predict(dpgvmat)
    assert glgmod.score(dpgvmat, bvmat) == 1.0

########################################
###### Variance calculation tests ######
########################################


########################################
######## Selection limit tests #########
########################################
def test_usl(glgmod, dpgvmat, mu, beta, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = beta.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1) + mu.flatten()
    b_usl = glgmod.usl(dpgvmat)

    assert numpy.all(a_usl == b_usl)

def test_usl(glgmod, dpgvmat, mu, beta, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = beta.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_lsl = uhaplo.sum(1) + mu.flatten()
    b_lsl = glgmod.lsl(dpgvmat)

    assert numpy.all(a_lsl == b_lsl)

### File I/O tests ###
def test_to_from_hdf5(glgmod, shared_datadir):
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "glgmod.hdf5")

    glgmod1 = GenericLinearGenomicModel.from_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod2 = GenericLinearGenomicModel.from_hdf5(
        shared_datadir / "glgmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(glgmod.mu == glgmod1.mu)
    assert numpy.all(glgmod.beta == glgmod1.beta)
    assert numpy.all(glgmod.trait == glgmod1.trait)
    assert glgmod.model_name == glgmod1.model_name
    assert glgmod.params == glgmod1.params

    assert numpy.all(glgmod.mu == glgmod2.mu)
    assert numpy.all(glgmod.beta == glgmod2.beta)
    assert numpy.all(glgmod.trait == glgmod2.trait)
    assert glgmod.model_name == glgmod2.model_name
    assert glgmod.params == glgmod2.params
