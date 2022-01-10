import numpy
import pytest
import os.path

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybropt.model.gmod.AdditiveLinearGenomicModel import is_AdditiveLinearGenomicModel
from pybropt.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from pybropt.model.gmod.AdditiveLinearGenomicModel import cond_check_is_AdditiveLinearGenomicModel

from pybropt.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybropt.popgen.bvmat.BreedingValueMatrix import is_BreedingValueMatrix
from pybropt.core.util import is_ndarray

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
def mat_u():
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
    yield "test_glgmod"

@pytest.fixture
def params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def glgmod(mat_beta, mat_u, mat_trait, model_name, params):
    yield AdditiveLinearGenomicModel(
        beta = mat_beta,
        u = mat_u,
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
##################### Breeding values ######################
############################################################
@pytest.fixture
def mat_intercept(dpgmat, mat_beta):
    n = dpgmat.ntaxa
    q = mat_beta.shape[0]
    a = numpy.ones((n,1), dtype = "float64")
    yield a

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(AdditiveLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "__init__")

def test_copy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "__copy__")

def test_deepcopy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "__deepcopy__")

def test_fit_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "fit_numpy")

def test_fit_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "fit")

def test_predict_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "predict_numpy")

def test_predict_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "predict")

def test_score_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "score_numpy")

def test_score_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "score")

def test_gebv_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "gebv_numpy")

def test_gebv_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "gebv")

def test_var_G_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_G_numpy")

def test_var_G_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_G")

def test_var_A_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_A_numpy")

def test_var_A_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_A")

def test_var_a_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_a_numpy")

def test_var_a_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "var_a")

def test_bulmer_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "bulmer_numpy")

def test_bulmer_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "bulmer")

def test_usl_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "usl_numpy")

def test_usl_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "usl")

def test_lsl_numpy_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "lsl_numpy")

def test_lsl_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "lsl")

def test_from_hdf5_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "from_hdf5")

def test_to_hdf5_is_concrete():
    generic_assert_concrete_method(AdditiveLinearGenomicModel, "to_hdf5")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_u_fget(glgmod, mat_u):
    assert numpy.all(glgmod.u == mat_u)

def test_beta_fget(glgmod, mat_beta):
    assert numpy.all(glgmod.beta == mat_beta)

def test_trait_fget(glgmod, mat_trait):
    assert numpy.all(glgmod.trait == mat_trait)

def test_model_name_fget(glgmod, model_name):
    assert glgmod.model_name == model_name

def test_params_fget(glgmod, params):
    assert glgmod.params == params

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########################################
########### Prediction tests ###########
########################################
def test_fit_numpy(glgmod):
    with pytest.raises(AttributeError):
        glgmod.fit_numpy(None, None, None)

def test_fit(glgmod):
    with pytest.raises(AttributeError):
        glgmod.fit(None, None, None)

def test_predict_numpy(glgmod, mat_intercept, mat_beta, mat_int8, mat_u):
    geno = mat_int8.sum(0)
    a = glgmod.predict_numpy(mat_intercept, geno)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u)
    assert numpy.all(a == b)
    assert is_ndarray(a)

def test_predict(glgmod, mat_intercept, mat_beta, mat_int8, mat_u, dpgmat):
    geno = mat_int8.sum(0)
    a = glgmod.predict(mat_intercept, dpgmat)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u)
    b = (b - b.mean(0)) / b.std(0)
    assert numpy.all(a == b)
    assert is_BreedingValueMatrix(a)

def test_score_numpy(glgmod, mat_intercept, mat_int8):
    geno = mat_int8.sum(0)
    y_true = glgmod.predict_numpy(mat_intercept, geno)
    out = glgmod.score_numpy(y_true, mat_intercept, geno)
    assert is_ndarray(out)
    assert numpy.all(out == 1.0)
    assert len(out) == glgmod.ntrait

def test_score(glgmod, mat_intercept, dpgmat):
    y_true = glgmod.predict(mat_intercept, dpgmat)
    out = glgmod.score(y_true, mat_intercept, dpgmat)
    assert is_ndarray(out)
    assert numpy.all(out == 1.0)

########################################
###### Variance calculation tests ######
########################################


########################################
######## Selection limit tests #########
########################################
def test_usl(glgmod, dpgmat, mat_u, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1)
    b_usl = glgmod.usl(dpgmat)

    assert numpy.all(a_usl == b_usl)

def test_usl(glgmod, dpgmat, mat_u, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)    pass

    a_lsl = uhaplo.sum(1)
    b_lsl = glgmod.lsl(dpgmat)

    assert numpy.all(a_lsl == b_lsl)

### File I/O tests ###
def test_to_from_hdf5(glgmod, shared_datadir):
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod.to_hdf5(shared_datadir / "glgmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "glgmod.hdf5")

    glgmod1 = AdditiveLinearGenomicModel.from_hdf5(shared_datadir / "glgmod.hdf5")
    glgmod2 = AdditiveLinearGenomicModel.from_hdf5(
        shared_datadir / "glgmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(glgmod.beta == glgmod1.beta)
    assert numpy.all(glgmod.u == glgmod1.u)
    assert numpy.all(glgmod.trait == glgmod1.trait)
    assert glgmod.model_name == glgmod1.model_name
    assert glgmod.params == glgmod1.params

    assert numpy.all(glgmod.beta == glgmod2.beta)
    assert numpy.all(glgmod.u == glgmod2.u)
    assert numpy.all(glgmod.trait == glgmod2.trait)
    assert glgmod.model_name == glgmod2.model_name
    assert glgmod.params == glgmod2.params

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_AdditiveLinearGenomicModel_is_concrete():
    generic_assert_concrete_function(is_AdditiveLinearGenomicModel)

def test_check_is_AdditiveLinearGenomicModel_is_concrete():
    generic_assert_concrete_function(check_is_AdditiveLinearGenomicModel)

def test_cond_check_is_AdditiveLinearGenomicModel_is_concrete():
    generic_assert_concrete_function(cond_check_is_AdditiveLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_AdditiveLinearGenomicModel(glgmod):
    assert is_AdditiveLinearGenomicModel(glgmod)

def test_check_is_AdditiveLinearGenomicModel(glgmod):
    with not_raises(TypeError):
        check_is_AdditiveLinearGenomicModel(glgmod, "glgmod")
    with pytest.raises(TypeError):
        check_is_AdditiveLinearGenomicModel(None, "glgmod")
