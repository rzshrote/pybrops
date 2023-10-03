import numpy
import pytest
import os.path

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import check_is_DenseAdditiveLinearGenomicModel

from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

################################ Test fixtures #################################

################# General shape parameters #################
@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def ntaxa():
    yield 100

@pytest.fixture
def nmarker():
    yield 100

@pytest.fixture
def nphase():
    yield 2

@pytest.fixture
def nchrom():
    yield 10

@pytest.fixture
def ngroup():
    yield 5

###################### Genomic model #######################
@pytest.fixture
def algmod_beta(ntrait):
    out = numpy.random.uniform(0, 10, (1,ntrait))
    yield out

@pytest.fixture
def algmod_u_misc():
    yield None

@pytest.fixture
def algmod_u_a(nmarker, ntrait):
    out = numpy.random.normal(size = (nmarker, ntrait))
    yield out

@pytest.fixture
def algmod_trait(ntrait):
    out = numpy.array(["Trait"+str(i+1) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def algmod_model_name():
    yield "test_algmod"

@pytest.fixture
def algmod_params():
    yield {"a" : 0, "b" : 1}

@pytest.fixture
def algmod(
        algmod_beta, 
        algmod_u_misc, 
        algmod_u_a, 
        algmod_trait, 
        algmod_model_name, 
        algmod_params
    ):
    yield DenseAdditiveLinearGenomicModel(
        beta       = algmod_beta,
        u_misc     = algmod_u_misc,
        u_a        = algmod_u_a,
        trait      = algmod_trait,
        model_name = algmod_model_name,
        params     = algmod_params
    )

######################## Genotypes #########################
@pytest.fixture
def pgmat_mat(nphase, ntaxa, nmarker):
    out = numpy.random.randint(0,2,(nphase,ntaxa,nmarker))
    out = out.astype("int8")
    yield out

@pytest.fixture
def pgmat_chrgrp(nchrom, nmarker):
    out = numpy.random.randint(1, nchrom+1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_phypos(nmarker):
    out = numpy.random.randint(1, 2**32-1, nmarker)
    out.sort()
    yield out

@pytest.fixture
def pgmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i+1).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pgmat_taxa_grp(ngroup, ntaxa):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pgmat(
        pgmat_mat, 
        pgmat_chrgrp, 
        pgmat_phypos, 
        pgmat_taxa, 
        pgmat_taxa_grp
    ):
    yield DensePhasedGenotypeMatrix(
        mat = pgmat_mat,
        vrnt_chrgrp = pgmat_chrgrp,
        vrnt_phypos = pgmat_phypos,
        taxa = pgmat_taxa,
        taxa_grp = pgmat_taxa_grp
    )

##################### Breeding values ######################
@pytest.fixture
def mat_intercept(pgmat, algmod_beta):
    n = pgmat.ntaxa
    q = algmod_beta.shape[0]
    a = numpy.ones((n,1), dtype = "float64")
    yield a

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(DenseAdditiveLinearGenomicModel)

############################ Test Class Properties #############################

def test_u_misc_fget(algmod, algmod_u_misc):
    assert numpy.all(algmod.u_misc == algmod_u_misc)

def test_u_a_fget(algmod, algmod_u_a):
    assert numpy.all(algmod.u_a == algmod_u_a)

def test_beta_fget(algmod, algmod_beta):
    assert numpy.all(algmod.beta == algmod_beta)

def test_trait_fget(algmod, algmod_trait):
    assert numpy.all(algmod.trait == algmod_trait)

def test_algmod_model_name_fget(algmod, algmod_model_name):
    assert algmod.model_name == algmod_model_name

def test_algmod_params_fget(algmod, algmod_params):
    assert algmod.params == algmod_params

########################## Test Class Special Methods ##########################

### __init__
def test___init___is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__init__")

### __copy__
def test___copy___is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__copy__")

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "__deepcopy__")

############################# Test concrete methods ############################

########### Prediction tests ###########

### fit_numpy
def test_fit_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "fit_numpy")

def test_fit_numpy(algmod):
    with pytest.raises(AttributeError):
        algmod.fit_numpy(None, None, None)

### fit
def test_fit_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "fit")

def test_fit(algmod):
    with pytest.raises(AttributeError):
        algmod.fit(None, None, None)

### predict_numpy
def test_predict_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "predict_numpy")

def test_predict_numpy(algmod, mat_intercept, algmod_beta, pgmat_mat, algmod_u_a):
    geno = pgmat_mat.sum(0)
    a = algmod.predict_numpy(mat_intercept, geno)
    b = (mat_intercept @ algmod_beta) + (geno @ algmod_u_a)
    assert numpy.all(a == b)
    assert isinstance(a, numpy.ndarray)

### predict
def test_predict_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "predict")

# def test_predict(algmod, mat_intercept, algmod_beta, pgmat_mat, algmod_u_a, pgmat):
#     geno = pgmat_mat.sum(0)
#     a = algmod.predict(mat_intercept, pgmat)
#     b = (mat_intercept @ algmod_beta) + (geno @ algmod_u_a)
#     b = (b - b.mean(0)) / b.std(0)
#     assert numpy.all(a == b)
#     assert isinstance(a, BreedingValueMatrix)

### score_numpy
def test_score_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "score_numpy")

def test_score_numpy(algmod, mat_intercept, pgmat_mat):
    geno = pgmat_mat.sum(0)
    y_true = algmod.predict_numpy(mat_intercept, geno)
    out = algmod.score_numpy(y_true, mat_intercept, geno)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)
    assert len(out) == algmod.ntrait

### score
def test_score_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "score")

def test_score(algmod, mat_intercept, pgmat):
    y_true = algmod.predict(mat_intercept, pgmat)
    out = algmod.score(y_true, mat_intercept, pgmat)
    assert isinstance(out, numpy.ndarray)
    assert numpy.all(out == 1.0)

### gebv_numpy
def test_gebv_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "gebv_numpy")

### gebv
def test_gebv_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "gebv")

###### Variance calculation tests ######

### var_G_numpy
def test_var_G_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_G_numpy")

### var_G
def test_var_G_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_G")

### var_A_numpy
def test_var_A_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_A_numpy")

### var_A
def test_var_A_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_A")

### var_a_numpy
def test_var_a_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_a_numpy")

### var_a
def test_var_a_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "var_a")

### bulmer_numpy
def test_bulmer_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "bulmer_numpy")

### bulmer
def test_bulmer_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "bulmer")

######## Selection limit tests #########

### usl_numpy
def test_usl_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "usl_numpy")

### usl
def test_usl_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "usl")

def test_usl(algmod, pgmat, algmod_u_a, pgmat_mat):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = algmod_u_a.T[:,None,None,:] * pgmat_mat
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.max((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)
    a_usl = uhaplo.sum(1)
    b_usl = algmod.usl(pgmat)

    assert numpy.all(a_usl == b_usl)

### lsl_numpy
def test_lsl_numpy_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "lsl_numpy")

### lsl
def test_lsl_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "lsl")

def test_lsl(algmod, pgmat, algmod_u_a, pgmat_mat):
    # (p,t)' -> (t,p)
    # (t,p)[:,None,None,:] -> (t,1,1,p)
    # (t,1,1,p) * (m,n,p) -> (t,m,n,p)
    haplo = algmod_u_a.T[:,None,None,:] * pgmat_mat
    # (t,m,n,p).max[1,2] -> (t,p)
    # scalar * (t,p) -> (t,p)
    uhaplo = 2.0 * haplo.min((1,2))
    # (t,p).sum(1) -> (t,)
    # (t,1).flatten -> (t,)
    # (t,) + (t,) -> (t,)    pass

    a_lsl = uhaplo.sum(1)
    b_lsl = algmod.lsl(pgmat)

    assert numpy.all(a_lsl == b_lsl)

############ File I/O tests ############

def test_to_hdf5_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "to_hdf5")

def test_from_hdf5_is_concrete():
    assert_concrete_method(DenseAdditiveLinearGenomicModel, "from_hdf5")

def test_to_from_hdf5(algmod, shared_datadir):
    algmod.to_hdf5(shared_datadir / "algmod.hdf5")
    algmod.to_hdf5(shared_datadir / "algmod.hdf5", "prefix")

    # test whether file was created
    assert os.path.isfile(shared_datadir / "algmod.hdf5")

    algmod1 = DenseAdditiveLinearGenomicModel.from_hdf5(shared_datadir / "algmod.hdf5")
    algmod2 = DenseAdditiveLinearGenomicModel.from_hdf5(
        shared_datadir / "algmod.hdf5",
        "prefix"
    )

    # test whether data was loaded properly
    assert numpy.all(algmod.beta == algmod1.beta)
    assert numpy.all(algmod.u_misc == algmod1.u_misc)
    assert numpy.all(algmod.u_a == algmod1.u_a)
    assert numpy.all(algmod.trait == algmod1.trait)
    assert algmod.model_name == algmod1.model_name
    assert algmod.params == algmod1.params

    assert numpy.all(algmod.beta == algmod2.beta)
    assert numpy.all(algmod.u_misc == algmod2.u_misc)
    assert numpy.all(algmod.u_a == algmod2.u_a)
    assert numpy.all(algmod.trait == algmod2.trait)
    assert algmod.model_name == algmod2.model_name
    assert algmod.params == algmod2.params

######################### Test class utility functions #########################

### check_is_DenseAdditiveLinearGenomicModel
def test_check_is_DenseAdditiveLinearGenomicModel_is_concrete():
    assert_concrete_function(check_is_DenseAdditiveLinearGenomicModel)

def test_check_is_DenseAdditiveLinearGenomicModel(algmod):
    with not_raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(algmod, "algmod")
    with pytest.raises(TypeError):
        check_is_DenseAdditiveLinearGenomicModel(None, "algmod")
