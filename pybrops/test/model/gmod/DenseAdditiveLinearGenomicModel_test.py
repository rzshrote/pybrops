import numpy
import pytest
import os.path

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import is_DenseAdditiveLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import check_is_DenseAdditiveLinearGenomicModel

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import is_BreedingValueMatrix

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
    assert_docstring(DenseAdditiveLinearGenomicModel)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__init__")

def test_copy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__copy__")

def test_deepcopy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__deepcopy__")

def test_fit_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "fit_numpy")

def test_fit_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "fit")

def test_predict_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "predict_numpy")

def test_predict_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "predict")

def test_score_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "score_numpy")

def test_score_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "score")

def test_gebv_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "gebv_numpy")

def test_gebv_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "gebv")

def test_var_G_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_G_numpy")

def test_var_G_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_G")

def test_var_A_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_A_numpy")

def test_var_A_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_A")

def test_var_a_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_a_numpy")

def test_var_a_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_a")

def test_bulmer_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "bulmer_numpy")

def test_bulmer_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "bulmer")

def test_usl_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "usl_numpy")

def test_usl_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "usl")

def test_lsl_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "lsl_numpy")

def test_lsl_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "lsl")

def test_from_hdf5_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "from_hdf5")

def test_to_hdf5_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "to_hdf5")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_u_misc_fget(dalgmod, mat_u_misc):
    assert numpy.all(dalgmod.u_misc == mat_u_misc)

def test_u_a_fget(dalgmod, mat_u_a):
    assert numpy.all(dalgmod.u_a == mat_u_a)

def test_beta_fget(dalgmod, mat_beta):
    assert numpy.all(dalgmod.beta == mat_beta)

def test_trait_fget(dalgmod, mat_trait):
    assert numpy.all(dalgmod.trait == mat_trait)

def test_model_name_fget(dalgmod, model_name):
    assert dalgmod.model_name == model_name

def test_params_fget(dalgmod, params):
    assert dalgmod.params == params

################################################################################
###################### Test concrete method functionality ######################
################################################################################

########################################
########### Prediction tests ###########
########################################
def test_fit_numpy(dalgmod):
    with pytest.raises(AttributeError):
        dalgmod.fit_numpy(None, None, None)

def test_fit(dalgmod):
    with pytest.raises(AttributeError):
        dalgmod.fit(None, None, None)

def test_predict_numpy(dalgmod, mat_intercept, mat_beta, mat_int8, mat_u_a):
    geno = mat_int8.sum(0)
    a = dalgmod.predict_numpy(mat_intercept, geno)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u_a)
    assert numpy.all(a == b)
    assert isinstance(a, numpy.ndarray)

def test_predict(dalgmod, mat_intercept, mat_beta, mat_int8, mat_u_a, dpgmat):
    geno = mat_int8.sum(0)
    a = dalgmod.predict(mat_intercept, dpgmat)
    b = (mat_intercept @ mat_beta) + (geno @ mat_u_a)
    b = (b - b.mean(0)) / b.std(0)
    assert numpy.all(a == b)
    assert is_BreedingValueMatrix(a)

def test_score_numpy(dalgmod, mat_intercept, mat_int8):
    geno = mat_int8.sum(0)
    y_true = dalgmod.predict_numpy(mat_intercept, geno)
    out = dalgmod.score_numpy(y_true, mat_intercept, geno)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)
    assert len(out) == dalgmod.ntrait

def test_score(dalgmod, mat_intercept, dpgmat):
    y_true = dalgmod.predict(mat_intercept, dpgmat)
    out = dalgmod.score(y_true, mat_intercept, dpgmat)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)

########################################
###### Variance calculation tests ######
########################################


########################################
######## Selection limit tests #########
########################################
def test_usl(dalgmod, dpgmat, mat_u_a, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u_a.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1)
    b_usl = dalgmod.usl(dpgmat)

    assert numpy.all(a_usl == b_usl)

def test_usl(dalgmod, dpgmat, mat_u_a, mat_int8):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = mat_u_a.T[:,None,None,:] * mat_int8
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)    pass

    a_lsl = uhaplo.sum(1)
    b_lsl = dalgmod.lsl(dpgmat)

    assert numpy.all(a_lsl == b_lsl)

### File I/O tests ###
def test_to_from_hdf5(dalgmod, shared_datadir):
    dalgmod.to_hdf5(shared_datadir / "dalgmod.hdf5")
    dalgmod.to_hdf5(shared_datadir / "dalgmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "dalgmod.hdf5")

    dalgmod1 = DenseAdditiveLinearGenomicModel.from_hdf5(shared_datadir / "dalgmod.hdf5")
    dalgmod2 = DenseAdditiveLinearGenomicModel.from_hdf5(
        shared_datadir / "dalgmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(dalgmod.beta == dalgmod1.beta)
    assert numpy.all(dalgmod.u_misc == dalgmod1.u_misc)
    assert numpy.all(dalgmod.u_a == dalgmod1.u_a)
    assert numpy.all(dalgmod.trait == dalgmod1.trait)
    assert dalgmod.model_name == dalgmod1.model_name
    assert dalgmod.params == dalgmod1.params

    assert numpy.all(dalgmod.beta == dalgmod2.beta)
    assert numpy.all(dalgmod.u_misc == dalgmod2.u_misc)
    assert numpy.all(dalgmod.u_a == dalgmod2.u_a)
    assert numpy.all(dalgmod.trait == dalgmod2.trait)
    assert dalgmod.model_name == dalgmod2.model_name
    assert dalgmod.params == dalgmod2.params

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DenseAdditiveLinearGenomicModel_is_concrete():
    assert_concrete_function(is_DenseAdditiveLinearGenomicModel)

def test_check_is_DenseAdditiveLinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_DenseAdditiveLinearGenomicModel)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseAdditiveLinearGenomicModel(dalgmod):
    assert is_DenseAdditiveLinearGenomicModel(dalgmod)

def test_check_is_DenseAdditiveLinearGenomicModel(dalgmod):
    with not_raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(dalgmod, "dalgmod")
    with pytest.raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(None, "dalgmod")
